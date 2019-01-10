import pytest
from reaktoro import *


@pytest.fixture()
def create_minimum_system():
    """
    Build a system with only a aqueous phase
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O")

    system = ChemicalSystem(editor)

    return system


def verify_dict(chemical):
    dict1 = {'teste':chemical}
    dict1[chemical] = 'teste'

    dict2 = {chemical:'teste'}
    dict2['teste'] = chemical
    assert dict1['teste'] == dict2['teste'] and dict1[chemical] == dict2[chemical]


@pytest.mark.xfail
def test_chemical_state_with_empty_system():
    system = ChemicalSystem()
    state = ChemicalState(system)
    repr(state)


def test_dict_on_chemical_state(create_minimum_system):
    state = ChemicalState(create_minimum_system)
    verify_dict(state)


def test_dict_on_chemical_quantity(create_minimum_system):
    quantity = ChemicalQuantity(create_minimum_system)
    verify_dict(quantity)
