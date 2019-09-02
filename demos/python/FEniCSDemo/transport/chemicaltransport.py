# Python modules
from dolfin import *
import time as timer
import numpy
from reaktoro.PyReaktoro import Partition, EquilibriumSolver, EquilibriumOptions
from transport.transport import TransportSolver



class ChemicalDirichletBC:
    def __init__(self, function_space, state, mobility, boundary):
        self.fluid = mobility.fluid(state)
        self.values = Function(function_space)
        self.dirichlet = DirichletBC(function_space, 1, boundary)
        self.dirichlet.apply(self.values.vector())
        self.dofs = numpy.where(self.values.vector() == 1)[0]
        self.dirichlet = DirichletBC(function_space, self.values, boundary)

    def elementDirichletBC(self, ielement, ispecies, porosity):
        bval = self.fluid.elementAmountInSpecies(ielement, ispecies)
        phi = porosity.vector()[self.dofs]
        self.values.vector()[self.dofs] = phi*bval
        return self.dirichlet

    def speciesDirichletBC(self, ispecies, porosity):
        nval = self.fluid.speciesAmount(ispecies)
        phi = porosity.vector()[self.dofs]
        self.values.vector()[self.dofs] = phi*nval
        return self.dirichlet


class ChemicalTransportResult(object):
    class EquilibriumResult(object):
        def __init__(self):
            # The number of iterations for equilibrium calculations performed at every degree of freedom
            self.iterations = []

            # The elapsed seconds for equilibrium calculations performed at every degree of freedom
            self.seconds = []

    class KineticsResult(object):
        def __init__(self):
            # The number of time-steps for kinetic calculations performed at every degree of freedom
            self.timesteps = []

            # The elapsed seconds for kinetic calculations performed at every degree of freedom
            self.seconds = []

    def __init__(self):
        # The flag that indicates if the chemical transport calculations succeeded
        self.succeeded = False

        # The total time in seconds of the chemical transport calculation
        self.seconds = 0

        # The total time in seconds spent on equilibrium calculations
        self.seconds_equilibrium = 0

        # The total time in seconds spent on finite element assembly operations
        self.seconds_assemble = 0

        # The total time in seconds spent on linear systems
        self.seconds_linear_systems = 0

        # The result of the equilibrium calculations
        self.equilibrium = ChemicalTransportResult.EquilibriumResult()

        # The result of the kinetic calculations
        self.kinetics = ChemicalTransportResult.KineticsResult()


class ChemicalTransportSolver(object):

    def __init__(self, field):
        self.system = field.system()
        self.mobility = field.mobility()
        self.partition = Partition(self.system)
        self.num_species  = self.system.numSpecies()
        self.num_elements = self.system.numElements()
        self.num_fluid_phases = len(self.mobility.indicesFluidPhases())
        self.num_solid_phases = len(self.mobility.indicesSolidPhases())
        self.velocity = [Constant(0.0) for i in range(self.num_fluid_phases)]
        self.diffusion = [Constant(0.0) for i in range(self.num_fluid_phases)]
        self.source = [Constant(0.0) for i in range(self.num_fluid_phases)]
        self.boundary_conditions = []
        self.initialized = False

        # Initialize the function space used in the ChemicalField instance
        self.function_space = field.functionSpace()

        # Initialize the list of global `dof` indices in the current process
        self.dofs = self.function_space.dofmap().dofs()

        # Initialize the number of degrees-of-freedom in the current process
        self.num_dofs = len(self.dofs)

        # Initialize the ChemicalTransportResult instance
        self.result = ChemicalTransportResult()
        self.result.equilibrium.iterations = Function(self.function_space)
        self.result.equilibrium.seconds = Function(self.function_space)
        self.result.kinetics.timesteps = Function(self.function_space)
        self.result.kinetics.seconds = Function(self.function_space)

        self.result.equilibrium._iterations = numpy.empty(self.num_dofs)
        self.result.equilibrium._seconds = numpy.empty(self.num_dofs)
        self.result.kinetics._timesteps = numpy.empty(self.num_dofs)
        self.result.kinetics._seconds = numpy.empty(self.num_dofs)

        self.result.equilibrium.iterations.rename('EquilibriumIterationsPerDOF', 'EquilibriumIterationsPerDOF')
        self.result.equilibrium.seconds.rename('EquilibriumSecondsPerDOF', 'EquilibriumSecondsPerDOF')
        self.result.kinetics.timesteps.rename('KineticsTimeStepsPerDOF', 'KineticsTimeStepsPerDOF')
        self.result.kinetics.seconds.rename('KineticsSecondsPerDOF', 'KineticsSecondsPerDOF')


    def setVelocity(self, velocity):
        assert type(velocity) is list, \
            'Expecting a list of velocity fields for each fluid phase.'
        assert len(velocity) == self.num_fluid_phases, \
            'There are %d fluid phases, but only %d velocity fields were given.' \
                % (self.num_fluid_phases, len(velocity))
        self.velocity = velocity
        self.initialized = False


    def setDiffusion(self, diffusion):
        assert type(diffusion) is list, \
            'Expecting a list of diffusion coefficients for each fluid phase.'
        assert len(diffusion) == self.num_fluid_phases, \
            'There are %d fluid phases, but only %d diffusion coefficients were given.' \
                % (self.num_fluid_phases, len(diffusion))
        self.diffusion = diffusion
        self.initialized = False


    def setSource(self, source):
        self.source = source
        self.initialized = False


    def setPartition(self, partition):
        if type(partition) is str:
            self.partition.set(partition)
        else:
            self.partition = partition


    def addBoundaryCondition(self, state, boundary):
        self.boundary_conditions.append((state, boundary))


    def setTemperatures(self, temperatures):
        if type(temperatures) in [int, float]:
            temperatures = [temperatures]*self.num_dofs
        assert hasattr(temperatures, '__len__'), 'Expecting a list of temperatures.'
        assert len(temperatures) == self.num_dofs, 'Expecting a list of \
            temperatures of length %d, but given list has length %d.' \
                % (self.num_dofs, len(temperatures))
        for i in range(self.num_dofs):
            self.states[i].setTemperature(temperatures[i])


    def setPressures(self, pressures):
        if type(pressures) in [int, float]:
            pressures = [pressures]*self.num_dofs
        assert hasattr(pressures, '__len__'), 'Expecting a list of pressures.'
        assert len(pressures) == self.num_dofs, 'Expecting a list of \
            pressures of length %d, but given list has length %d.' \
                % (self.num_dofs, len(pressures))
        for i in range(self.num_dofs):
            self.states[i].setPressure(pressures[i])


    def initialize(self, field):

        self.initialized = True

        # Check if the user provided any boundary conditions
        if self.boundary_conditions is []:
            RuntimeError('Failed to initialize ChemicalTransportSolver. \
                No boundary conditions have been provided.')

        # Initialize the ChemicalDirichletBC instances
        self.bcs = [ChemicalDirichletBC(self.function_space, state, self.mobility, boundary) \
                    for (state, boundary) in self.boundary_conditions]

        # The auxiliary Function instance used for transport steps
        self.u = Function(self.function_space)

        # The auxiliary Function instance used for outputting
        self.output = Function(self.function_space)

        # Initialize the chemical equilibrium solver
        self.equilibrium = EquilibriumSolver(self.system)
        self.equilibrium.setPartition(self.partition)

        # Initialize the indices of the equilibrium and kinetic species
        self.ispecies_e  = self.partition.indicesEquilibriumSpecies()
        self.ispecies_k  = self.partition.indicesKineticSpecies()

        # Initialize the indices of the fluid and solid species
        self.ispecies_f  = self.mobility.indicesFluidSpeciesInEachFluidPhase()
        self.ispecies_s  = self.mobility.indicesSolidSpecies()

        # Initialize the indices of the equilibrium-fluid and equilibrium-solid species
        self.ispecies_ef = [sorted(set(self.ispecies_e) & set(indices)) for indices in self.ispecies_f]
        self.ispecies_es = sorted(set(self.ispecies_e) & set(self.ispecies_s))

        # Initialize the indices of the kinetic-fluid and kinetic-solid species
        self.ispecies_kf = [sorted(set(self.ispecies_k) & set(indices)) for indices in self.ispecies_f]
        self.ispecies_ks = sorted(set(self.ispecies_k) & set(self.ispecies_s))

        # Initialize the indices of fluid and solid phases
        self.iphases_f = self.mobility.indicesFluidPhases()
        self.iphases_s = self.mobility.indicesSolidPhases()

        # Initialize the number of fluid and solid phases
        self.num_fluid_phases = len(self.iphases_f)
        self.num_solid_phases = len(self.iphases_s)

        # Initialize the arrays of element amounts for each fluid phase in
        # the equilibrium-fluid partition for each degree-of-freedom in the function space
        self.bef = [numpy.zeros((self.num_elements, self.num_dofs)) \
                        for i in range(self.num_fluid_phases)]

        # Initialize the array of element amounts in the equilibrium
        # partition for each degree-of-freedom in the function space
        self.be = numpy.zeros((self.num_elements, self.num_dofs))

        # Get the dolfin Function's for the saturation fields of fluid phases
        self.saturations = field.saturations()

        # Get the dolfin Function for the porosity field
        self.porosity = field.porosity()

        # Define the dolfin forms for the pore velocity of each fluid phase
        self.pore_velocity = [self.velocity[i]/(self.porosity*self.saturations[i]) \
            for i in range(self.num_fluid_phases)]

        # Create transport solvers for each fluid phase
        self.transport = [TransportSolver() for i in range(self.num_fluid_phases)]

        # Set the pore velocities and diffusion coefficients for each transport solver
        for i in range(self.num_fluid_phases):
            self.transport[i].setVelocity(self.pore_velocity[i])
            self.transport[i].setDiffusion(self.diffusion[i])

        # Initialize the porosity and saturations fields
        field.update()

        # Initialize the element amounts in the fluid phases of the equilibrium-fluid partition
        for i in range(self.num_fluid_phases):
            field.elementAmountsInSpecies(self.ispecies_ef[i], self.bef[i])


    def elementDirichletBC(self, ielement, iphase):
        ispecies = self.ispecies_ef[iphase] # the indices of equilibrium species in phase `iphase`
        return [bc.elementDirichletBC(ielement, ispecies, self.porosity) for bc in self.bcs]


    def speciesDirichletBC(self, ispecies):
        return [bc.speciesDirichletBC(ispecies, self.porosity) for bc in self.bcs]


    def transportEquilibriumElementsInFluidPhases(self, field, dt):
        # Iterate over all elements in all fluid phases in the equilibrium-fluid
        # partition and tranport them
        for iphase in range(self.num_fluid_phases):
            
            # Iterate over all elements in the current fluid phase and tranport therm
            for ielement in range(self.num_elements):
                # The Dirichlet boundary conditions element `ielement` in fluid phase `iphase`
                bcs = self.elementDirichletBC(ielement, iphase)

                # Set the boundary conditions in the tranport solver w.r.t. fluid phase `iphase`
                self.transport[iphase].setBoundaryConditions(bcs)

                # Set initial condition for the transport equation of the current element
                self.u.vector()[:] = self.bef[iphase][ielement]

                # Transport the current element of the equilibrium-fluid partition
                self.transport[iphase].step(self.u, dt)

                # Extract the result to the array of element amounts `bef`
                self.bef[iphase][ielement][:] = self.u.vector().get_local()


    def transportKineticSpeciesInFluidPhases(self, field, dt):
        pass


    def equilibrate(self, field):

        # Compute the positive and negative part of the molar amounts of the elements in the equilibrium-fluid partition
        bef_positive = [numpy.maximum(x, 0) for x in self.bef]
        bef_negative = [numpy.minimum(x, 0) for x in self.bef]

        # Start timing the equilibrate step
        tbegin = timer.time()

        # Define auxiliary bindings
        states = field.states()
        iterations = self.result.equilibrium._iterations
        seconds = self.result.equilibrium._seconds

        # Calculate the element amounts in the equilibrium partition.
        # Compute the contribution from the equilibrium-solid partition
        field.elementAmountsInSpecies(self.ispecies_es, self.be)

        # Compute the contribution from the equilibrium-fluid partition
        for bef in bef_positive:
            self.be += bef

        options = EquilibriumOptions()
#         options.epsilon = 1e-30
#         options.hessian = EquilibriumHessian.Exact
#         options.optimum.tolerance = 1e-8
#         options.optimum.output.active = True

        self.equilibrium.setOptions(options)

        # Compute the equilibrium states
        for k in range(self.num_dofs):
            
            # Perform the equilibrium calculation at current degree of freedom
            result = self.equilibrium.solve(states[k], states[k].temperature(),states[k].pressure(), self.be[:, k])

            # Check if the equilibrium calculation was successful
            if not result.optimum.succeeded:
                b = [(elem.name(), amount) for elem, amount in zip(self.system.elements(), self.be[:, k])]
                raise RuntimeError("Failed to calculate equilibrium state at dof (%d), under " \
                "temperature %f K, pressure %f Pa, and element molar amounts %s." % \
                (k, states[k].temperature(), states[k].pressure(), str(b)))

            # Store the statistics of the equilibrium calculation
            iterations[k] = result.optimum.iterations
            seconds[k] = result.optimum.time

        # Update the element amounts in the fluid phases of the equilibrium-fluid partition
        for i in range(self.num_fluid_phases):
            field.elementAmountsInSpecies(self.ispecies_ef[i], self.bef[i])
            self.bef[i] += bef_negative[i]

        # Extract the calculation statistics from the auxiliary ndarray members
        self.result.equilibrium.iterations.vector()[:] = iterations
        self.result.equilibrium.seconds.vector()[:] = seconds

        # Total time spent on performing equilibrium calculations
        self.result.seconds_equilibrium = timer.time() - tbegin

    def react(self, field, dt):
        pass

    def step(self, field, dt):
        # Check if the chemical transport solver has been initialized
        if not self.initialized:
            self.initialize(field)

        # Start timing the transport step calculation
        tbegin = timer.time()

        # Transport the elements in the equilibrium-fluid partition
        self.transportEquilibriumElementsInFluidPhases(field, dt)

        # Transport the species in the kinetic-fluid partition
        self.transportKineticSpeciesInFluidPhases(field, dt)

        # Equilibrate the equilibrium-solid and equilibrium-fluid partitions
        self.equilibrate(field)

        # React the equilibrium and kinetic partitions over a time of `dt`
        self.react(field, dt)

        # Update the porosity and fluid phases saturation fields
        field.update()

        # Compute the time elapsed for the chemical transport step
        self.result.time = timer.time() - tbegin

    def elementAmountInPhase(self, element, phase):
        ielement = self.system.indexElement(element)
        iphase = self.system.indexPhase(phase)
        out = Function(self.function_space)
        out.vector()[:] = self.bef[iphase][ielement]
        out.rename(element, element)
        return out