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
#include <functional>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/Outputter.hpp>

namespace Reaktoro {

/// A type that describes the evaluation result of a non-linear residual function.
struct NonlinearResidual
{
    /// The residual vector of the non-linear residual function evaluated at `x`.
    Vector val;

    /// The Jacobian matrix of the non-linear residual function evaluated at `x`.
    Matrix jacobian;

    /// The boolean flag that indicates if the non-linear function evaluation succeeded.
    /// This should be set to false whenever the evaluation of the non-linear function
    /// fails. For example, when it produces `inf` or `nan` values, or when this function
    /// relies on more complicated ones that can fail.
    bool succeeded = true;
};

/// A type that describes the functional signature of a non-linear residual function.
/// @param x The vector of variables
/// @return The residual of the non-linear function evaluated at `x`
using NonlinearFunction = std::function<NonlinearResidual(VectorConstRef x)>;

/// A type that describes the non-linear problem.
struct NonlinearProblem
{
    /// The non-linear residual function.
    NonlinearFunction f;

    /// The number of unknowns in the non-linear problem.
    Index n;

    /// The number of non-linear functions.
    Index m;

    /// The left-hand side matrix of the linear inequality constraint `A*x = b`.
    /// This can be left empty if no linear inequality constraints exist.
    Matrix A;

    /// The right-hand side vector of the linear equality constraint `A*x = b`.
    /// This can be left empty if no linear inequality constraints exist.
    Vector b;
};

/// A type that describes the options for the output of a non-linear problem calculation.
struct NonlinearOutput : OutputterOptions
{
    /// The prefix for the unknown variables `x`.
    std::string xprefix = "x";

    /// The prefix for the residual functions.
    std::string fprefix = "f";

    /// The names of the unknown variables `x`.
    /// Numbers will be used if not properly set (e.g., `x[0]`, `x[1]`)
    std::vector<std::string> xnames;

    /// The names of the non-linear functions.
    /// Numbers will be used if not properly set (e.g., `f[0]`, `f[1]`)
    std::vector<std::string> fnames;

    /// Activate or deactivate output.
    auto operator=(bool activate) -> NonlinearOutput&;
};

/// A type that describes the options for the solution of a non-linear problem.
struct NonlinearOptions
{
    /// The tolerance for the residual of the non-linear function.
    double tolerance = 1.0e-6;

    /// The tolerance for the variation in variables x.
    /// Set this to a value greater than zero to stop the calculation
    /// whenever `max(abs(dx)) < tolerancex`, where `dx` is the current step
    /// of the unknown variables.
    double tolerancex = 0.0;

    /// The maximum number of iterations in the solution of the non-linear problem.
    unsigned max_iterations = 100;

    /// The boundary to the fraction parameter.
    double tau = 0.9999;

    /// The Armijo parameter used in the backtracking line search algorithm.
    double armijo = 1.0e-4;

    /// The flag that indicates if backtracking line search algorithm is to be used.
    bool linesearch = true;

    /// The options for the output of the non-linear problem calculation.
    NonlinearOutput output;
};

/// A type that describes the result of a non-linear problem calculation.
struct NonlinearResult
{
    /// The flag that indicates if the calculation converged.
    bool succeeded = false;

    /// The number of iterations in the calculation.
    unsigned iterations = 0;

    /// The number of evaluations of the residual non-linear function in the calculation.
    unsigned num_function_evals = 0;

    /// The final residual error of the calculation.
    double error = 0;

    /// The wall time spent for the calculation (in units of s).
    double time = 0;

    /// The wall time spent for all function evaluations (in units of s).
    double time_function_evals = 0;

    /// The wall time spent for all linear system solutions (in units of s).
    double time_linear_systems = 0;
};

/// A type that implements the Newton algorithm for solving non-linear problems.
class NonlinearSolver
{
public:
    /// Construct a default NonlinearSolver instance.
    NonlinearSolver();

    /// Construct a copy of an NonlinearSolver instance.
    NonlinearSolver(const NonlinearSolver& other);

    /// Destroy this NonlinearSolver instance.
    virtual ~NonlinearSolver();

    /// Assign an NonlinearSolver instance to this.
    auto operator=(NonlinearSolver other) -> NonlinearSolver&;

    /// Solve a non-linear problem.
    /// @param problem The definition of the non-linear problem.
    /// @param x[in,out] The initial guess and the final solution of the calculation.
    auto solve(const NonlinearProblem& problem, VectorRef x) -> NonlinearResult;

    /// Solve a non-linear problem with given options.
    /// @param problem The definition of the non-linear problem.
    /// @param x[in,out] The initial guess and the final solution of the calculation.
    /// @param options The options for the calculation.
    auto solve(const NonlinearProblem& problem, VectorRef x, const NonlinearOptions& options) -> NonlinearResult;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
