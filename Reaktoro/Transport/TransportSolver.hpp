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
#include <memory>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
//class ChemicalProperties;
//class ChemicalState;
//class ChemicalSystem;

//class BoundaryState
//{
//public:
//    BoundaryState(const ChemicalState& state);
//
//private:
//};
//
//class BoundaryFlux
//{
//public:
//    BoundaryFlux(const ChemicalState& state);
//
//private:
//};


class ChemicalField
{
public:
    using Iterator = std::vector<ChemicalState>::iterator;

    using ConstIterator = std::vector<ChemicalState>::const_iterator;

    ChemicalField(Index size, const ChemicalSystem& system);

    ChemicalField(Index size, const ChemicalState& state);

    auto size() const -> Index { return m_size; }

    auto begin() const -> ConstIterator { return m_states.cbegin(); }

    auto begin() -> Iterator { return m_states.begin(); }

    auto end() const -> ConstIterator { return m_states.cend(); }

    auto end() -> Iterator { return m_states.end(); }

    auto operator[](Index index) const -> const ChemicalState& { return m_states[index]; }

    auto operator[](Index index) -> ChemicalState& { return m_states[index]; }

    auto set(const ChemicalState& state) -> void;

    auto temperature(VectorRef values) -> void;

    auto pressure(VectorRef values) -> void;

    auto elementAmounts(VectorRef values) -> void;

    auto output(std::string filename, StringList quantities) -> void;

private:
    /// The number of degrees of freedom in the chemical field.
    Index m_size;

//    Vector temperatures;
//
//    Vector pressures;
//
//    /// The matrix of amounts for every element (
//    Matrix element_amounts;

    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem m_system;

    /// The chemical states in the chemical field
    std::vector<ChemicalState> m_states;

    /// The chemical states in the chemical field
    std::vector<ChemicalProperties> m_properties;
};

/// A class that defines a Tridiagonal Matrix used on TransportSolver.
/// it stores data in a Eigen::VectorXd like, M = {a[0][0], a[0][1], a[0][2], 
///                                                a[1][0], a[1][1], a[1][2],
///                                                a[2][0], a[2][1], a[2][2]}
class TridiagonalMatrix
{
public:
    TridiagonalMatrix() : TridiagonalMatrix(0) {}

    TridiagonalMatrix(Index size) : m_size(size), m_data(size * 3) {}

    auto size() const -> Index { return m_size; }

    auto data() -> VectorRef { return m_data; }

    auto data() const -> VectorConstRef { return m_data; }

    auto row(Index index) -> VectorRef { return m_data.segment(3 * index, 3); }

    auto row(Index index) const -> VectorConstRef { return m_data.segment(3 * index, 3); }

    auto a() -> VectorStridedRef { return Vector::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    auto a() const -> VectorConstRef { return Vector::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    auto b() -> VectorStridedRef { return Vector::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    auto b() const -> VectorConstRef { return Vector::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    auto c() -> VectorStridedRef { return Vector::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    auto c() const -> VectorConstRef { return Vector::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    auto resize(Index size) -> void;

    /// Factorize the tridiagonal matrix to help with linear system A x = b.
    auto factorize() -> void;

    /// Solve a linear system Ax = b with LU decomposition.
    auto solve(VectorRef x, VectorConstRef b) const -> void;

    /// Solve a linear system Ax = b with LU decomposition, using x as the unknown and it's.
    /// old values as the vector b.
    auto solve(VectorRef x) const -> void;

    operator Matrix() const;

private:
    /// The size of the tridiagonal matrix
    Index m_size;

    /// The coefficients
    Vector m_data;
};

/// A class that defines the mesh for TransportSolver.
class Mesh
{
public:
    Mesh();

    Mesh(Index num_cells, double xl = 0.0, double xr = 1.0);

    auto setDiscretization(Index num_cells, double xl = 0.0, double xr = 1.0) -> void;

    auto numCells() const -> Index { return m_num_cells; }

    auto xl() const -> double { return m_xl; }

    auto xr() const -> double { return m_xr; }

    auto dx() const -> double { return m_dx; }

    auto xcells() const -> VectorConstRef { return m_xcells; }

private:
    /// The number of cells in the discretization.
    Index m_num_cells = 10;

    /// The x-coordinate of the left boundary (in m).
    double m_xl = 0.0;

    /// The x-coordinate of the right boundary (in m).
    double m_xr = 1.0;

    /// The length of the cells (in m).
    double m_dx = 0.1;

    /// The x-coordinate of the center of the cells.
    Vector m_xcells;
};

/// A class for solving advection-diffusion problem.
/// Eq: du/dt + v*du/dx = D*d�u/dx�
///     u - amount
///     v - velocity
///     D - diffusion coefficient
class TransportSolver
{
public:
    /// Construct a default TransportSolver instance.
    TransportSolver();

    /// Set the mesh for the numerical solution of the transport problem.
    auto setMesh(const Mesh& mesh) -> void { mesh_ = mesh; }

    /// Set the velocity for the transport problem.
    /// @param val The velocity (in m/s)
    auto setVelocity(double val) -> void { velocity = val; }

    /// Set the diffusion coefficient for the transport problem.
    /// @param val The diffusion coefficient (in m^2/s)
    auto setDiffusionCoeff(double val) -> void { diffusion = val; }

    /// Set the value of the variable on the boundary.
    /// @param val The boundary value for the variable (same unit considered for u).
    auto setBoundaryValue(double val) -> void { ul = val; };

    /// Set the time step for the numerical solution of the transport problem.
    auto setTimeStep(double val) -> void { dt = val; }

    /// Return the mesh.
    auto mesh() const -> const Mesh& { return mesh_; }

    /// Initialize the transport solver before method @ref step is executed.
    /// Setup coefficient matrix of the diffusion problem and factorize.
    auto initialize() -> void;

    /// Step the transport solver.
    /// This method solve one step of the transport solver equation, using an explicit approach for
    /// advection and total implicit for diffusion. The amount resulted from the advection it is 
    /// passed to diffusion problem as a "source".
    /// @param[in,out] u The solution vector
    /// @param q The source rates vector ([same unit considered for u]/m)
    auto step(VectorRef u, VectorConstRef q) -> void;

    /// Step the transport solver.
    /// @param[in,out] u The solution vector
    auto step(VectorRef u) -> void;

private:
    /// The mesh describing the discretization of the domain.
    Mesh mesh_;

    /// The time step used to solve the transport problem (in s).
    double dt = 0.0;

    /// The velocity in the transport problem (in m/s).
    double velocity = 0.0;

    /// The diffusion coefficient in the transport problem (in m^2/s).
    double diffusion = 0.0;

    /// The value of the variable on the left boundary.
    double ul;

    /// The coefficient matrix from the discretized transport equation
    TridiagonalMatrix A;

    /// The flux limiters at each cell.
    Vector phi;

    /// The previous state of the variables.
    Vector u0;
};

/// Use this class for solving reactive transport problems.
class ReactiveTransportSolver
{
public:
    /// Construct a default ReactiveTransportSolver instance.
    ReactiveTransportSolver(const ChemicalSystem& system);

    auto setMesh(const Mesh& mesh) -> void;

    auto setVelocity(double val) -> void;

    auto setDiffusionCoeff(double val) -> void;

    auto setBoundaryState(const ChemicalState& state) -> void;

    auto setTimeStep(double val) -> void;

    auto system() const -> const ChemicalSystem& { return system_; }

    auto output() -> ChemicalOutput;

    auto initialize() -> void;

    auto step(ChemicalField& field) -> void;

private:
    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem system_;

    /// The solver for solving the transport equations
    TransportSolver transportsolver;

    /// The solver for solving the equilibrium equations
    EquilibriumSolver equilibriumsolver;

    /// The list of chemical output objects
    std::vector<ChemicalOutput> outputs;

    /// The amounts of fluid elements on the boundary.
    Vector bbc;

    /// The amounts of a fluid element on each cell of the mesh.
    Matrix bf;

    /// The amounts of a solid element on each cell of the mesh.
    Matrix bs;

    /// The amounts of an element on each cell of the mesh.
    Matrix b;

    /// The current number of steps in the solution of the reactive transport equations.
    Index steps = 0;
};

} // namespace Reaktoro
