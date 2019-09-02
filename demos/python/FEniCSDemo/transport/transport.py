from dolfin import *


class _TransportSolver(object):
    def __init__(self):

        # Initialize default values for the velocity, diffusion and source parameters
        self.velocity = Constant(0.0)
        self.diffusion = Constant(0.0)
        self.source = Constant(0.0)

        # Initialize the list of DirichletBC instances
        self.bcs = []

        # The time step as a dolfin.Constant instance to avoid compilation when its value changes
        self.dt = Constant(0.0)

        # The flag that indicates if the solver has been initialized before evolving the solution
        self.initialized = False

    def setVelocity(self, velocity):
        self.velocity = velocity
        self.initialized = False

    def setDiffusion(self, diffusion):
        self.diffusion = diffusion
        self.initialized = False

    def setSource(self, source):
        self.source = source
        self.initialized = False

    def setBoundaryConditions(self, bcs):
        self.bcs = bcs if hasattr(bcs, '__len__') else [bcs]

    def initialize(self, u):
        self.initialized = True

        velocity = self.velocity
        diffusion = self.diffusion
        source = self.source
        dt = self.dt

        V = u.function_space()
        mesh = V.mesh()

        u0 = u

        u = TrialFunction(V)
        v = TestFunction(V)

        # Mid-point solution
        u_mid = 0.5*(u0 + u)

        # Residual
        r = u - u0 + dt*(dot(velocity, grad(u_mid)) - \
            diffusion*div(grad(u_mid)) - source)

        # Galerkin variational problem
        F = v*(u - u0)*dx + dt*(v*div(u_mid*velocity)*dx + \
            diffusion*dot(grad(v), grad(u_mid))*dx - source*v*dx)

        # Add SUPG stabilisation terms
        h = 2*Circumradius(mesh)
        vnorm = sqrt(dot(velocity, velocity))
        tau = h/(2.0*vnorm)

        F += tau*dot(velocity, grad(v))*r*dx

        # Create bilinear and linear forms
        self.a = lhs(F)
        self.L = rhs(F)

        self.A = PETScMatrix()
        self.b = PETScVector()

        self.solver = LUSolver()
        self.solver.parameters["symmetric"] = True

    def step(self, u, dt):
        if not self.initialized:
            self.initialize(u)
        self.dt.assign(dt)
        assemble(self.a, tensor=self.A)
        assemble(self.L, tensor=self.b)
        for bc in self.bcs:
            bc.apply(self.A)
            bc.apply(self.b)
        self.solver.set_operator(self.A)
        self.solver.solve(u.vector(), self.b)


class TransportSolver(object):
    def __init__(self):
        self.pimpl = _TransportSolver()

    def setVelocity(self, velocity):
        self.pimpl.setVelocity(velocity)

    def setDiffusion(self, diffusion):
        self.pimpl.setDiffusion(diffusion)

    def setSource(self, source):
        self.pimpl.setSource(source)

    def setBoundaryConditions(self, bc):
        self.pimpl.setBoundaryConditions(bc)

    def step(self, u, dt):
        self.pimpl.step(u, dt)
