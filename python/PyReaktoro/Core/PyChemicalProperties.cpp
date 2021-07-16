// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

void exportChemicalProperties(py::module& m)
{
    auto update1 = static_cast<void (ChemicalProperties::*)(double, double)>(&ChemicalProperties::update);
    auto update2 = static_cast<void (ChemicalProperties::*)(VectorConstRef)>(&ChemicalProperties::update);
    auto update3 = static_cast<void (ChemicalProperties::*)(double, double, VectorConstRef)>(&ChemicalProperties::update);
    auto update4 = static_cast<void (ChemicalProperties::*)(double, double, VectorConstRef, const ThermoModelResult&, const ChemicalModelResult&)>(&ChemicalProperties::update);

    py::class_<ChemicalProperties>(m, "ChemicalProperties")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("update", update1)
        .def("update", update2)
        .def("update", update3)
        .def("update", update4)
        .def("temperature", &ChemicalProperties::temperature)
        .def("pressure", &ChemicalProperties::pressure)
        .def("composition", &ChemicalProperties::composition)
        .def("thermoModelResult", &ChemicalProperties::thermoModelResult, py::return_value_policy::reference_internal)
        .def("chemicalModelResult", &ChemicalProperties::chemicalModelResult, py::return_value_policy::reference_internal)
        .def("moleFractions", &ChemicalProperties::moleFractions)
        .def("lnActivityCoefficients", &ChemicalProperties::lnActivityCoefficients, py::return_value_policy::reference_internal)
        .def("lnActivityConstants", &ChemicalProperties::lnActivityConstants, py::return_value_policy::reference_internal)
        .def("lnActivities", &ChemicalProperties::lnActivities, py::return_value_policy::reference_internal)
        .def("partialMolarVolumes", &ChemicalProperties::partialMolarVolumes, py::return_value_policy::reference_internal)
        .def("chemicalPotentials", &ChemicalProperties::chemicalPotentials)
        .def("standardPartialMolarGibbsEnergies", &ChemicalProperties::standardPartialMolarGibbsEnergies, py::return_value_policy::reference_internal)
        .def("standardPartialMolarEnthalpies", &ChemicalProperties::standardPartialMolarEnthalpies, py::return_value_policy::reference_internal)
        .def("standardPartialMolarVolumes", &ChemicalProperties::standardPartialMolarVolumes, py::return_value_policy::reference_internal)
        .def("standardPartialMolarEntropies", &ChemicalProperties::standardPartialMolarEntropies)
        .def("standardPartialMolarInternalEnergies", &ChemicalProperties::standardPartialMolarInternalEnergies)
        .def("standardPartialMolarHelmholtzEnergies", &ChemicalProperties::standardPartialMolarHelmholtzEnergies)
        .def("standardPartialMolarHeatCapacitiesConstP", &ChemicalProperties::standardPartialMolarHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("standardPartialMolarHeatCapacitiesConstV", &ChemicalProperties::standardPartialMolarHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        .def("phaseMolarGibbsEnergies", &ChemicalProperties::phaseMolarGibbsEnergies)
        .def("phaseMolarEnthalpies", &ChemicalProperties::phaseMolarEnthalpies)
        .def("phaseMolarVolumes", &ChemicalProperties::phaseMolarVolumes)
        .def("phaseMolarEntropies", &ChemicalProperties::phaseMolarEntropies)
        .def("phaseMolarInternalEnergies", &ChemicalProperties::phaseMolarInternalEnergies)
        .def("phaseMolarHelmholtzEnergies", &ChemicalProperties::phaseMolarHelmholtzEnergies)
        .def("phaseMolarHeatCapacitiesConstP", &ChemicalProperties::phaseMolarHeatCapacitiesConstP)
        .def("phaseMolarHeatCapacitiesConstV", &ChemicalProperties::phaseMolarHeatCapacitiesConstV)
        .def("phaseSpecificGibbsEnergies", &ChemicalProperties::phaseSpecificGibbsEnergies)
        .def("phaseSpecificEnthalpies", &ChemicalProperties::phaseSpecificEnthalpies)
        .def("phaseSpecificVolumes", &ChemicalProperties::phaseSpecificVolumes)
        .def("phaseSpecificEntropies", &ChemicalProperties::phaseSpecificEntropies)
        .def("phaseSpecificInternalEnergies", &ChemicalProperties::phaseSpecificInternalEnergies)
        .def("phaseSpecificHelmholtzEnergies", &ChemicalProperties::phaseSpecificHelmholtzEnergies)
        .def("phaseSpecificHeatCapacitiesConstP", &ChemicalProperties::phaseSpecificHeatCapacitiesConstP)
        .def("phaseSpecificHeatCapacitiesConstV", &ChemicalProperties::phaseSpecificHeatCapacitiesConstV)
        .def("phaseDensities", &ChemicalProperties::phaseDensities)
        .def("phaseMasses", &ChemicalProperties::phaseMasses)
        .def("phaseAmounts", &ChemicalProperties::phaseAmounts)
        .def("phaseVolumes", &ChemicalProperties::phaseVolumes)
        .def("phaseInternalEnergies", &ChemicalProperties::phaseInternalEnergies)
        .def("phaseEnthalpies", &ChemicalProperties::phaseEnthalpies)
        .def("volume", &ChemicalProperties::volume)
        .def("internalEnergy", &ChemicalProperties::internalEnergy)
        .def("enthalpy", &ChemicalProperties::enthalpy)
        .def("subvolume", &ChemicalProperties::subvolume)
        .def("fluidVolume", &ChemicalProperties::fluidVolume)
        .def("solidVolume", &ChemicalProperties::solidVolume)
        ;
}

} // namespace Reaktoro
