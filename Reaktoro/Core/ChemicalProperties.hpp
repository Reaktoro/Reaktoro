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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Thermodynamics/Models/ChemicalModel.hpp>
#include <Reaktoro/Thermodynamics/Models/ThermoModel.hpp>

namespace Reaktoro {

/// A class for querying thermodynamic and chemical properties of a chemical system.
class ChemicalProperties
{
public:
    /// Construct a default ChemicalProperties instance.
    ChemicalProperties();

    /// Construct a ChemicalProperties instance with given ChemicalSystem.
    ChemicalProperties(const ChemicalSystem& system);

    /// Update the thermodynamic properties of the chemical system.
    /// @param T The temperature in the system (in units of K)
    /// @param P The pressure in the system (in units of Pa)
    auto update(double T, double P) -> void;

    /// Update the chemical properties of the chemical system.
    /// @param n The amounts of the species in the system (in units of mol)
    auto update(VectorConstRef n) -> void;

    /// Update the thermodynamic and chemical properties of the chemical system.
    /// @param T The temperature in the system (in units of K)
    /// @param P The pressure in the system (in units of Pa)
    /// @param n The amounts of the species in the system (in units of mol)
    auto update(double T, double P, VectorConstRef n) -> void;

    /// Update the thermodynamic and chemical properties of the chemical system.
    /// @param T The temperature in the system (in units of K)
    /// @param P The pressure in the system (in units of Pa)
    /// @param n The amounts of the species in the system (in units of mol)
    /// @param tres The result of the ThermoModel function of the chemical system.
    /// @param cres The result of the ChemicalModel function of the chemical system.
    auto update(double T, double P, VectorConstRef n, const ThermoModelResult& tres, const ChemicalModelResult& cres) -> void;

    /// Return the temperature of the system (in units of K).
    auto temperature() const -> Temperature;

    /// Return the pressure of the system (in units of Pa).
    auto pressure() const -> Pressure;

    /// Return the molar amounts of the species (in units of mol).
    auto composition() const -> Composition;

    /// Return the result of the PhaseThermoModel function of each phase.
    auto thermoModelResult() const -> const ThermoModelResult&;

    /// Return the result of the PhaseChemicalModel function of each phase.
    auto chemicalModelResult() const -> const ChemicalModelResult&;

    /// Return the mole fractions of the species.
    auto moleFractions() const -> ChemicalVector;

    /// Return the ln activity coefficients of the species.
    auto lnActivityCoefficients() const -> ChemicalVectorConstRef;

    /// Return the ln activity constants of the species.
    auto lnActivityConstants() const -> ThermoVectorConstRef;

    /// Return the ln activities of the species.
    auto lnActivities() const -> ChemicalVectorConstRef;

    /// Return the partial molar volume of the species in that phase (in units of m3/mol).
    auto partialMolarVolumes() const -> ChemicalVectorConstRef;

    /// Return the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials() const -> ChemicalVector;

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> ThermoVectorConstRef;

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> ThermoVectorConstRef;

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> ThermoVectorConstRef;

    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropies() const -> ThermoVector;

    /// Return the standard partial molar internal energies of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergies() const -> ThermoVector;

    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector;

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVectorConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVectorConstRef;

    /// Return the molar Gibbs energies of the phases (in units of J/mol).
    auto phaseMolarGibbsEnergies() const -> ChemicalVector;

    /// Return the molar enthalpies of the phases (in units of J/mol).
    auto phaseMolarEnthalpies() const -> ChemicalVector;

    /// Return the molar volumes of the phases (in units of m3/mol).
    auto phaseMolarVolumes() const -> ChemicalVector;

    /// Return the molar entropies of the phases (in units of J/(mol*K)).
    auto phaseMolarEntropies() const -> ChemicalVector;

    /// Return the molar internal energies of the phases (in units of J/mol).
    auto phaseMolarInternalEnergies() const -> ChemicalVector;

    /// Return the molar Helmholtz energies of the phases (in units of J/mol).
    auto phaseMolarHelmholtzEnergies() const -> ChemicalVector;

    /// Return the molar isobaric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstP() const -> ChemicalVector;

    /// Return the molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstV() const -> ChemicalVector;

    /// Return the specific Gibbs energies of the phases (in units of J/kg).
    auto phaseSpecificGibbsEnergies() const -> ChemicalVector;

    /// Return the specific enthalpies of the phases (in units of J/kg).
    auto phaseSpecificEnthalpies() const -> ChemicalVector;

    /// Return the specific volumes of the phases (in units of m3/kg).
    auto phaseSpecificVolumes() const -> ChemicalVector;

    /// Return the specific entropies of the phases (in units of J/(kg*K)).
    auto phaseSpecificEntropies() const -> ChemicalVector;

    /// Return the specific internal energies of the phases (in units of J/kg).
    auto phaseSpecificInternalEnergies() const -> ChemicalVector;

    /// Return the specific Helmholtz energies of the phases (in units of J/kg).
    auto phaseSpecificHelmholtzEnergies() const -> ChemicalVector;

    /// Return the specific isobaric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector;

    /// Return the specific isochoric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector;

    /// Return the densities of the phases (in units of kg/m3).
    auto phaseDensities() const -> ChemicalVector;

    /// Return the masses of the phases (in units of kg).
    auto phaseMasses() const -> ChemicalVector;

    /// Return the molar amounts of the phases (in units of mol).
    auto phaseAmounts() const -> ChemicalVector;

    /// Return the volumes of the phases (in units of m3).
    auto phaseVolumes() const -> ChemicalVector;

    /// Return the internal energies of the phases (in units of J).
    auto phaseInternalEnergies() const -> ChemicalVector;

    /// Return the enthalpies of the phases (in units of J).
    auto phaseEnthalpies() const -> ChemicalVector;

    /// Return the volume of the system (in units of m3).
    auto volume() const -> ChemicalScalar;

    /// Return the internal energy of the system (in units of J).
    auto internalEnergy() const -> ChemicalScalar;

    /// Return the enthalpy of the system (in units of J).
    auto enthalpy() const -> ChemicalScalar;

    /// Return the total volume occupied by given phases (in units of m3).
    /// @param iphases The indices of the phases.
    auto subvolume(const Indices& iphases) const -> ChemicalScalar;

    /// Return the total fluid volume of the system (in units of m3).
    /// The fluid volume is defined as the sum of volumes of all fluid phases.
    auto fluidVolume() const -> ChemicalScalar;

    /// Return the total solid volume of the system (in units of m3).
    /// The solid volume is defined as the sum of volumes of all solid phases.
    auto solidVolume() const -> ChemicalScalar;

private:
    /// The chemical system
    ChemicalSystem system;

    /// The number of species in the system
    Index num_species;

    /// The number of phases in the system
    Index num_phases;

    /// The temperature of the system (in units of K)
    Temperature T;

    /// The pressure of the system (in units of Pa)
    Pressure P;

    /// The amounts of the species in the system (in units of mol).
    Vector n;

    /// The mole fractions of the species in the system (in units of mol/mol).
    ChemicalVector x;

    /// The results of the evaluation of the PhaseThermoModel functions of each phase.
    ThermoModelResult tres;

    /// The results of the evaluation of the PhaseChemicalModel functions of each phase.
    ChemicalModelResult cres;
};

} // namespace Reaktoro
