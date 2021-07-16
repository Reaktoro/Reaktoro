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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>

namespace Reaktoro {

/// A type used to describe the result of a smart equilibrium calculation.
struct SmartEquilibriumResult
{
    /// The boolean flag that indicates if smart equilibrium calculation was used.
    bool succeeded = false;
};

/// A type used to describe the result of an inverse equilibrium calculation.
struct InverseEquilibriumResult
{
    /// The number of forward equilibrium calculations in each inverse iteration.
    std::vector<Index> fcep_iterations_per_icep_iteration;

    /// The unknown variables *x* in each inverse iteration.
    std::vector<Vector> x_per_icep_iteration;

    /// The residual function *F* in each inverse iteration.
    std::vector<Vector> F_per_icep_iteration;

    /// The residual error *E* in each inverse iteration.
    std::vector<double> E_per_icep_iteration;
};

/// A type used to describe the result of an equilibrium calculation
/// @see ChemicalState
struct EquilibriumResult
{
    /// The result of the optimisation calculation
    OptimumResult optimum;

    /// The boolean flag that indicates if smart equilibrium calculation was used.
    SmartEquilibriumResult smart;

    /// The result details of an inverse equilibrium calculation.
    InverseEquilibriumResult inverse;

    /// Apply an addition assignment to this instance
    auto operator+=(const EquilibriumResult& other) -> EquilibriumResult&;
};

} // namespace Reaktoro
