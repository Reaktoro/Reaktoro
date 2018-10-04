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
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

// Forward declarations
class GaseousMixture;

/// Return an equation of state for a gaseous phase based on the ideal model.
/// @param mixture The gaseous mixture
/// @return The equation of state function for the gaseous phase
/// @see GaseousMixture, GaseousChemicalModel
auto gaseousChemicalModelIdeal(const GaseousMixture& mixture) -> PhaseChemicalModel;

} // namespace Reaktoro
