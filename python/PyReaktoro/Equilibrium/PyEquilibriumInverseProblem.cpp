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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>

namespace Reaktoro {

void exportEquilibriumInverseProblem(py::module& m)
{
    auto setTemperature1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::setTemperature);
    auto setTemperature2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::setTemperature);

    auto setPressure1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::setPressure);
    auto setPressure2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::setPressure);

    auto fixSpeciesAmount1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesAmount);
    auto fixSpeciesAmount2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesAmount);

    auto fixSpeciesMass1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesMass);
    auto fixSpeciesMass2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesMass);

    auto fixSpeciesActivity1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double)>(&EquilibriumInverseProblem::fixSpeciesActivity);
    auto fixSpeciesActivity2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesActivity);
    auto fixSpeciesActivity3 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesActivity);

    auto fixSpeciesFugacity1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string)>(&EquilibriumInverseProblem::fixSpeciesFugacity);
    auto fixSpeciesFugacity2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(std::string, double, std::string, std::string)>(&EquilibriumInverseProblem::fixSpeciesFugacity);

    auto pH1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::pH);
    auto pH2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::pH);
    auto pH3 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string, std::string)>(&EquilibriumInverseProblem::pH);

    auto pE1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double)>(&EquilibriumInverseProblem::pE);
    auto pE2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::pE);

    auto Eh1 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string)>(&EquilibriumInverseProblem::Eh);
    auto Eh2 = static_cast<EquilibriumInverseProblem&(EquilibriumInverseProblem::*)(double, std::string, std::string)>(&EquilibriumInverseProblem::Eh);

    py::class_<EquilibriumInverseProblem>(m, "EquilibriumInverseProblem")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const EquilibriumInverseProblem&>())
        .def("setPartition", &EquilibriumInverseProblem::setPartition, py::return_value_policy::reference_internal)
        .def("setTemperature", setTemperature1, py::return_value_policy::reference_internal)
        .def("setTemperature", setTemperature2, py::return_value_policy::reference_internal)
        .def("setPressure", setPressure1, py::return_value_policy::reference_internal)
        .def("setPressure", setPressure2, py::return_value_policy::reference_internal)
        .def("setElementInitialAmounts", &EquilibriumInverseProblem::setElementInitialAmounts, py::return_value_policy::reference_internal)
        .def("add", &EquilibriumInverseProblem::add, py::return_value_policy::reference_internal)
        .def("fixSpeciesAmount", fixSpeciesAmount1, py::return_value_policy::reference_internal)
        .def("fixSpeciesAmount", fixSpeciesAmount2, py::return_value_policy::reference_internal)
        .def("fixSpeciesMass", fixSpeciesMass1, py::return_value_policy::reference_internal)
        .def("fixSpeciesMass", fixSpeciesMass2, py::return_value_policy::reference_internal)
        .def("fixSpeciesChemicalPotential", &EquilibriumInverseProblem::fixSpeciesChemicalPotential, py::return_value_policy::reference_internal)
        .def("fixSpeciesActivity", fixSpeciesActivity1, py::return_value_policy::reference_internal)
        .def("fixSpeciesActivity", fixSpeciesActivity2, py::return_value_policy::reference_internal)
        .def("fixSpeciesActivity", fixSpeciesActivity3, py::return_value_policy::reference_internal)
        .def("fixSpeciesFugacity", fixSpeciesFugacity1, py::return_value_policy::reference_internal)
        .def("fixSpeciesFugacity", fixSpeciesFugacity2, py::return_value_policy::reference_internal)
        .def("fixPhaseAmount", &EquilibriumInverseProblem::fixPhaseAmount, py::return_value_policy::reference_internal)
        .def("fixPhaseMass", &EquilibriumInverseProblem::fixPhaseMass, py::return_value_policy::reference_internal)
        .def("fixPhaseVolume", &EquilibriumInverseProblem::fixPhaseVolume, py::return_value_policy::reference_internal)
        .def("fixPhaseSetVolume", &EquilibriumInverseProblem::fixPhaseSetVolume, py::return_value_policy::reference_internal)
        .def("pH", pH1, py::return_value_policy::reference_internal)
        .def("pH", pH2, py::return_value_policy::reference_internal)
        .def("pH", pH3, py::return_value_policy::reference_internal)
        .def("pE", pE1, py::return_value_policy::reference_internal)
        .def("pE", pE2, py::return_value_policy::reference_internal)
        .def("Eh", Eh1, py::return_value_policy::reference_internal)
        .def("Eh", Eh2, py::return_value_policy::reference_internal)
	    .def("alkalinity", &EquilibriumInverseProblem::alkalinity, py::return_value_policy::reference_internal)
	    .def("fixVolume", &EquilibriumInverseProblem::fixVolume, py::return_value_policy::reference_internal)
	    .def("fixInternalEnergy", &EquilibriumInverseProblem::fixInternalEnergy, py::return_value_policy::reference_internal)
	    .def("fixEnthalpy", &EquilibriumInverseProblem::fixEnthalpy, py::return_value_policy::reference_internal)
        .def("unknownTemperature", &EquilibriumInverseProblem::unknownTemperature)
        .def("unknownPressure", &EquilibriumInverseProblem::unknownPressure)
        .def("unknownAmountOf", &EquilibriumInverseProblem::unknownAmountOf)
        .def("unknownAmountOfEither", &EquilibriumInverseProblem::unknownAmountOfEither)
        .def("system", &EquilibriumInverseProblem::system, py::return_value_policy::reference_internal)
        .def("partition", &EquilibriumInverseProblem::partition, py::return_value_policy::reference_internal)
        .def("temperature", &EquilibriumInverseProblem::temperature)
        .def("pressure", &EquilibriumInverseProblem::pressure)
        .def("numConstraints", &EquilibriumInverseProblem::numConstraints)
        .def("numUnknowns", &EquilibriumInverseProblem::numUnknowns)
        .def("numTitrants", &EquilibriumInverseProblem::numTitrants)
        .def("formulaMatrixTitrants", &EquilibriumInverseProblem::formulaMatrixTitrants)
        .def("unknownsCoefficientMatrix", &EquilibriumInverseProblem::unknownsCoefficientMatrix)
        .def("elementInitialAmounts", &EquilibriumInverseProblem::elementInitialAmounts)
        .def("titrantInitialAmounts", &EquilibriumInverseProblem::titrantInitialAmounts)
        ;
}

} // namespace Reaktoro
