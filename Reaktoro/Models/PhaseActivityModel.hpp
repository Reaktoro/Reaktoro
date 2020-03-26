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
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

// Forward declarations
class Species;
class ThermoEngine;

/// The base type for all thermodynamic activity models for phases.
/// @see PhaseStandardThermoModel
class PhaseActivityModel
{
public:
    /// Create the activity model of the phase.
    /// @param engine The thermodynamic engine from which additional data can be fetched.
    /// @param species The species that compose the phase.
    virtual auto create(const ThermoEngine& engine, const std::vector<Species>& species) -> PhaseActivityModelFn = 0;
};

} // namespace Reaktoro
