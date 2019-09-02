import numpy
from dolfin import *

class Saturation(object):

    def __init__(self, field):
        Function.__init__(self, field.functionSpace())
        self.field = field
        self.mobility = field.mobility()
        self.ifluid_phases = self.mobility.indicesFluidPhases()
        self.num_fluid_phases = len(self.ifluid_phases)
        self.dofmap = field.functionSpace().dofmap()
        self.num_dofs = len(self.dofmap.dofs())
        self.saturation = [Function(field.functionSpace()) for i in self.ifluid_phases]
        self.values = numpy.zeros((self.num_fluid_phases, self.num_dofs))


    def update(self):
        states = self.field.states()
        for k in xrange(self.num_dofs):
            properties = states[k].properties()
            v = properties.phaseVolumes().val
            fluid_volume = sum([v[i] for i in self.ifluid_phases])
            self.values[:, k] = [v[i]/fluid_volume for i in self.ifluid_phases]
        for i in xrange(self.num_fluid_phases):
            self.saturation[i].vector()[:] = self.values[i]


    def __getitem__(self, ifluid_phase):
        return self.saturation[ifluid_phase]

