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

#include "EquilibriumInverseSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Optimization/NonlinearSolver.hpp>

namespace Reaktoro {

struct EquilibriumInverseSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// The sensitivity of the equilibrium state
    EquilibriumSensitivity sensitivity;

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system), solver(system)
    {
        // Set the default partition as all species are in equilibrium
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        // Set the partition of the chemical system
        partition = partition_;

        // Set the partition of the equilibrium solver
        solver.setPartition(partition);
    }

    /// Solve an inverse equilibrium problem
    auto solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
    {
        // The accumulated equilibrium result of this inverse problem calculation
        EquilibriumResult result;

        // Define auxiliary variables from the inverse problem definition
        const Index Nt = problem.numTitrants();
        const Index Nc = problem.numConstraints();
        const Index Nx = problem.numUnknowns();

        const Indices ies = partition.indicesEquilibriumSpecies();
        const Indices iee = partition.indicesEquilibriumElements();
        const Index Ne = partition.numEquilibriumSpecies();
        const Index Ee = partition.numEquilibriumElements();
        const Index Nu = 2 + Ee;

        // The coefficient matrix C used in u = u0 + C*x
        const Matrix C = problem.unknownsCoefficientMatrix();

        // Check the number of columns in C correspond to the number of constraints.
        if(Nx != Nc)
            RuntimeError("Could not solve the inverse chemical equilibrium problem",
                "The number of unknowns and constraints need to be equal.");

        // The initial temperature and pressure for the calculation
        const double T0 = problem.temperature();
        const double P0 = problem.pressure();

        // The initial element amounts for the calculation
        const Vector b0 = problem.elementInitialAmounts();
        const Vector be0 = rows(b0, iee);

        // The vector u with initial values
        Vector u0(2 + Ee);
        u0 << T0, P0, be0;

        Vector u(2 + Ee);

        // Set the temperature and pressure of the chemical state
        state.setTemperature(T0);
        state.setPressure(P0);

        // Define auxiliary instances to avoid memory reallocation
        ChemicalProperties properties;
        ResidualEquilibriumConstraints res;
        NonlinearResidual nonlinear_residual;

        // Auxiliary references to the non-linear residual data
        auto& F = nonlinear_residual.val;
        auto& J = nonlinear_residual.jacobian;

        // Set the options and partition in the equilibrium solver
        solver.setOptions(options);
        solver.setPartition(partition);

        // Define the non-linear problem with inequality constraints
        NonlinearProblem nonlinear_problem;

        // Set the linear inequality constraints of the titrant molar amounts
        nonlinear_problem.n = Nc;
        nonlinear_problem.m = Nc;

        // Set the non-linear function of the non-linear problem
        nonlinear_problem.f = [&](VectorConstRef x) mutable
        {
            // Update vector u
            u = u0 + C * x;

            // Extract T, P, be from u
            const auto T = u[0];
            const auto P = u[1];
            const auto be = u.tail(Ee);

            if(T <= 0.0 || P <= 0.0)
            {
                nonlinear_residual.succeeded = false;
                return nonlinear_residual;
            }

            // Solve the equilibrium problem with update `be`
            auto fcep_res = solver.solve(state, T, P, be);

            result += fcep_res;
            result.inverse.fcep_iterations_per_icep_iteration.push_back(fcep_res.optimum.iterations);

            // Update the chemical properties
            properties = solver.properties();

            // Check if the equilibrium calculation converged
            if(!result.optimum.succeeded)
            {
                // If not, solve using cold start
                state.setSpeciesAmounts(0.0);
                auto fcep_res = solver.solve(state, T, P, be);

                result += fcep_res;
                result.inverse.fcep_iterations_per_icep_iteration.back() += fcep_res.optimum.iterations;
            }

            // Check if the function evaluation was successful
            nonlinear_residual.succeeded = result.optimum.succeeded;

            // Update the sensitivity of the equilibrium state
            sensitivity = solver.sensitivity();

            // Alias for the T, P, b sensitivities
            const auto& dndT = sensitivity.dndT;
            const auto& dndP = sensitivity.dndP;
            const auto& dndb = sensitivity.dndb;

            Matrix dndu = zeros(Ne, Nu);
            dndu.col(0) = dndT;
            dndu.col(1) = dndP;
            dndu.rightCols(Ee) = dndb;

            // Calculate the residuals of the equilibrium constraints
            res = problem.residualEquilibriumConstraints(x, state, properties);

            // Calculate the residual vector `F` and its Jacobian `J`
            F = res.val;
            J = res.ddx;

            // Add T and P derivatives contributions to the first Nt columns of the Jacobian matrix
            J += res.ddu * C;
            J += res.ddn * dndu * C;

            nonlinear_residual.succeeded = true;

            result.inverse.x_per_icep_iteration.push_back(x);
            result.inverse.F_per_icep_iteration.push_back(F);
            result.inverse.E_per_icep_iteration.push_back(norminf(F));

            return nonlinear_residual;
        };

        // Initialize the initial guess of the titrant amounts
        // Vector x = problem.titrantInitialAmounts();
        Vector x = problem.initialGuessUnknowns();

        // Replace zeros in x by small molar amounts
        x = (x.array() > 0.0).select(x, 1e-6);

        // Solve the non-linear problem with inequality constraints
        NonlinearSolver nonlinear_solver;
        nonlinear_solver.solve(nonlinear_problem, x, options.nonlinear);

        return result;
    }
};

EquilibriumInverseSolver::EquilibriumInverseSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumInverseSolver::EquilibriumInverseSolver(const EquilibriumInverseSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumInverseSolver::~EquilibriumInverseSolver()
{}

auto EquilibriumInverseSolver::operator=(EquilibriumInverseSolver other) -> EquilibriumInverseSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumInverseSolver::setOptions(const EquilibriumOptions& options) -> void
{
    pimpl->options = options;
}

auto EquilibriumInverseSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumInverseSolver::solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult
{
    return pimpl->solve(state, problem);
}

auto EquilibriumInverseSolver::sensitivity() -> EquilibriumSensitivity
{
    return pimpl->solver.sensitivity();
}

} // namespace Reaktoro
