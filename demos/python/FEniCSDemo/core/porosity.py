import numpy
from dolfin import *

class Porosity(Function):

    def __init__(self, field):
        Function.__init__(self, field.functionSpace())
        self.field = field
        self.mobility = field.mobility()
        self.isolid_phases = self.mobility.indicesSolidPhases()
        self.dofmap = field.functionSpace().dofmap()
        self.num_dofs = len(self.dofmap.dofs())
        self.values = numpy.empty(self.num_dofs)

    def update(self):
        states = self.field.states()
        for k in xrange(self.num_dofs):
            properties = states[k].properties()
            v = properties.phaseVolumes().val
            solid_volume = sum([v[i] for i in self.isolid_phases])
            self.values[k] = 1.0 - solid_volume
        self.vector()[:] = self.values
