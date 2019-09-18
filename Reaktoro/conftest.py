import numpy as np
import pytest

from reaktoro import ChemicalEditor, ChemicalSystem, Partition


@pytest.fixture
def chemical_system():
    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--".split())
    editor.addGaseousPhase("H2O(g) CO2(g)".split())
    editor.addMineralPhase("Graphite")
    return ChemicalSystem(editor)


@pytest.fixture
def partition_with_inert_gaseous_phase(chemical_system):
    partition = Partition(chemical_system)
    partition.setInertPhases(['Gaseous'])

    return partition


@pytest.fixture
def chemical_properties(chemical_system):
    # A sensible value for temperature (in K)
    T = 300

    # A sensible value for pressure (in Pa)
    P = 1e5

    # A sensible array of species amounts
    n = np.array([55, 1e-7, 1e-7, 0.1, 0.5, 0.01, 1.0, 0.001, 1.0])

    return chemical_system.properties(T, P, n)

