'''
Created on Jul 22, 2015

@author: allan
'''

from dolfin import *
import time as stimer
from core.mobility import Mobility
from core.chemicalfield import ChemicalField
from transport.chemicaltransport import ChemicalTransportSolver
from reaktoro.PyReaktoro import Database, ChemicalEditor, ChemicalSystem,\
    EquilibriumProblem, ChemicalState, equilibrate
from utils import paraview


class InletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

minute = 60
hour = 60*minute
day = 24*hour

# tend = 6e+4
tend = 12000
dt = 300
t = 0.0

# The velocity (in units of m/s)
v = Constant((1.15740741e-5, 0.0)) # equivalent to 1 m/day

# The diffusion coefficient (in units of m2/s)
D = Constant(1.0e-9)

# Initialise the mesh
xl, xr = 0.0, 1.0
yb, yt = 0.0, 0.005
ncells = 100

mesh = RectangleMesh(Point(xl, yb), Point(xr, yt), ncells, 1)

CG = FunctionSpace(mesh, 'CG', 1)

# Initialise the database
database = Database(r'demos/supcrt98.xml')

# Initialise the chemical editor
editor = ChemicalEditor(database)
# editor.addAqueousPhase('H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- O2(aq) H2(aq)')
editor.addAqueousPhase('H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3--')
editor.addMineralPhase('Calcite')
editor.addMineralPhase('Dolomite')
editor.addMineralPhase('Quartz')

# Initialise the chemical system
system = ChemicalSystem(editor)

# Initialise the initial chemical state condition
problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(100, 'celsius')
problem_ic.setPressure(300, 'bar')
problem_ic.add('H2O',   1.00, 'kg')
problem_ic.add('NaCl',  1.0, 'mol')
problem_ic.add('MgCl2', 1e-6, 'mol')
problem_ic.add('CaCl2', 1e-6, 'mol')
problem_ic.add('CO2',   0.01, 'mol')
problem_ic.add('SiO2',  10, 'mol')
problem_ic.add('CaCO3', 10, 'mol')
problem_ic.add('O2',   1e-8, 'mol')

problem_ic.add('H2O', 1, 'kg')
problem_ic.add('NaCl', 1, 'mol')
problem_ic.add('Calcite', 10, 'mol')
problem_ic.add('Quartz', 10, 'mol')

state_ic = ChemicalState(system)

equilibrate(state_ic, problem_ic)

print(state_ic)

state_ic.scalePhaseVolume('Aqueous', 0.50, "m3")
state_ic.scalePhaseVolume('Calcite', 0.01, "m3")
state_ic.scalePhaseVolume('Quartz',  0.49, "m3")


# Initialise the boundary chemical state condition
problem_bc = EquilibriumProblem(system)
problem_bc.setTemperature(100, 'celsius')
problem_bc.setPressure(300, 'bar')
problem_bc.add('H2O',   1.00, 'kg')
problem_bc.add('NaCl',  1.0e-6, 'mol')
problem_bc.add('MgCl2', 0.05, 'mol')
problem_bc.add('CaCl2', 0.01, 'mol')
problem_bc.add('CO2',   0.75, 'mol')
problem_bc.add('SiO2',  1e-8, 'mol')
problem_ic.add('O2',   1e-8, 'mol')

state_bc = ChemicalState(system)

equilibrate(state_bc, problem_bc)


# Scale the volumes of the initial and boundary condition states to 1 m3
state_ic.scaleVolume(1.0)
state_bc.scaleVolume(1.0)

print('state_ic\n', state_ic)
print('state_bc\n', state_bc)

# Initialise the mobility info of the chemical system
mobility = Mobility(system)
mobility.setFluidPhases(['Aqueous'])

field = ChemicalField(system, mobility, CG)
field.fill(state_ic)

transport = ChemicalTransportSolver(field)
transport.addBoundaryCondition(state_bc, InletBoundary())
transport.setVelocity([v])
transport.setDiffusion([D])

out_species = ['Ca++', 'Mg++', 'Calcite', 'Dolomite', 'CO2(aq)', 'HCO3-', 'Cl-', 'H2O(l)']
out_elements = ['H', 'O', 'C', 'Ca', 'Mg', 'Na', 'Cl']

# Define the name of the result file
filename = r'demos/nosupg-result-100C-300bar-dt-%d-cells-%d.xdmf' % (dt, ncells)

# Create the output file
file = XDMFFile(mesh.mpi_comm(), filename)
begin = stimer.time()
while t <= tend:
    print('Time: ', t, '({:.2%})'.format(t/tend))

    # For each selected species, output its molar amounts
    for species in out_species:
        file.write ( field.speciesAmount(species),t)

    result = transport.result

    file.write ( result.equilibrium.iterations, t)
    file.write ( result.equilibrium.seconds, t)
    file.write ( field.porosity(), t)
    file.write ( field.volume(), t)

    # file.write ( field.ph(), t)

    # Perform one transport step from `t` to `t + dt`
    transport.step(field, dt)

    # For each selected element, output its molar amounts
    for element in out_elements:
        file.write(transport.elementAmountInPhase(element, 'Aqueous'), t)

    # Update the current time
    t += dt

paraview.xdmf(filename)

print('Statistics:')
print('Total simulation time = ', stimer.time() - begin)
print('Total time spent on transport calculations = ', result.time)
print('Total time spent on equilibrium calculations = ', result.seconds_equilibrium, '({:.2%})'.format(result.seconds_equilibrium/result.time))
