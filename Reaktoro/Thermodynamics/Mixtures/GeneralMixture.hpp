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

// C++ includes
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// A type used to describe the state of a mixture
struct MixtureState
{
    /// The temperature of the mixture (in units of K)
    real T;

    /// The pressure of the mixture (in units of Pa)
    real P;

    /// The mole fractions of the species in the mixture and their partial derivatives
    VectorXr x;
};

/// Compare two MixtureState instances for equality
auto operator==(const MixtureState& l, const MixtureState& r) -> bool;

/// Provide a base of implementation for the mixture classes.
/// @ingroup Mixtures
class GeneralMixture
{
public:
    /// Construct a default GeneralMixture instance
    GeneralMixture();

    /// Construct a GeneralMixture instance with given species
    /// @param species The names of the species in the mixture
    explicit GeneralMixture(const std::vector<Species>& species);

    /// Destroy the instance
    virtual ~GeneralMixture();

    /// Set the name of the mixture.
    auto setName(std::string name) -> void;

    /// Return the number of species in the mixture
    auto numSpecies() const -> unsigned;

    /// Return the name of the mixture.
    auto name() const -> std::string;

    /// Return the species that compose the mixture
    /// @return The species that compose the mixture
    auto species() const -> const std::vector<Species>&;

    /// Return a species in the mixture
    /// @param index The index of the species
    /// @return The species with given index
    auto species(const Index& index) const -> const Species&;

    /// Return the index of a species in the mixture
    /// @param name The name of the species in the mixture
    /// @return The index of the species if found, or the number of species otherwise
    auto indexSpecies(const std::string& name) const -> Index;

    /// Return the index of the first species in the mixture with any of the given names.
    /// @param names The tentative names of the species in the mixture.
    /// @return The index of the species if found, or the number of species otherwise.
    auto indexSpeciesAny(const std::vector<std::string>& names) const -> Index;

    /// Return the names of the species in the mixture
    auto namesSpecies() const -> std::vector<std::string>;

    /// Return the charges of the species in the mixture
    auto chargesSpecies() const -> VectorXr;

    /// Calculates the mole fractions of the species and their partial derivatives
    /// @param n The molar abundance of the species (in units of mol)
    /// @return The mole fractions and their partial derivatives
    auto moleFractions(VectorXrConstRef n) const -> VectorXr;

    /// Calculate the state of the mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(real T, real P, VectorXrConstRef n) const -> MixtureState;

private:
    /// The name of mixture
    std::string _name;

    /// The species in the mixture
    std::vector<Species> _species;
};

} // namespace Reaktoro
