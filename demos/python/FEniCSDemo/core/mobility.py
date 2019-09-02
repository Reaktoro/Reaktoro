from reaktoro.PyReaktoro import Connectivity

class Mobility(object):
    def __init__(self, system):
        self.system = system
        self.connectivity = Connectivity(self.system)
        self._ifluid_phases  = []
        self._isolid_phases  = range(system.numPhases())
        self._ifluid_species = []
        self._isolid_species = range(system.numSpecies())

    def setFluidPhases(self, fluid_phases):
        for phase_name in fluid_phases:
            iphase = self.system.indexPhaseWithError(phase_name)
            self._ifluid_phases.append(iphase)
            self._ifluid_species += self.connectivity.indicesSpeciesInPhase(iphase)
        prange = range(self.system.numPhases())
        srange = range(self.system.numSpecies())
        self._isolid_phases = set(prange) - set(self._ifluid_phases)
        self._isolid_species = set(srange) - set(self._ifluid_species)

    def setSolidPhases(self, solid_phases):
        fluid_phases = [phase.name() for phase in self.system.phases()]
        fluid_phases = set(fluid_phases) - set(solid_phases)
        self.setFluidPhases(fluid_phases)

    def indicesFluidPhases(self):
        return self._ifluid_phases

    def indicesFluidSpecies(self):
        return self._ifluid_species

    def indicesFluidSpeciesInEachFluidPhase(self):
        return [self.connectivity.indicesSpeciesInPhase(iphase) \
            for iphase in self.indicesFluidPhases()]

    def indicesSolidPhases(self):
        return self._isolid_phases

    def indicesSolidSpecies(self):
        return self._isolid_species

    def indicesSolidSpeciesInEachSolidPhase(self):
        return [self.connectivity.indicesSpeciesInPhase(iphase) \
            for iphase in self.indicesSolidPhases()]

    def fluid(self, state):
        res = state.clone()
        for i in self.indicesSolidPhases():
            res.scalePhaseVolume(i, 0.0)
        res.scaleVolume(1.0)
        return res

    def solid(self, state):
        res = state.clone()
        for i in self.indicesFluidPhases():
            res.scalePhaseVolume(i, 0.0)
        res.scaleVolume(1.0)
        return res

