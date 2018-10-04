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

namespace Reaktoro {

// Forward declarations
class ChemicalEditor;
class ChemicalPlot;
class ChemicalSystem;
class EquilibriumPath;
class EquilibriumProblem;
class KineticPath;
class MineralReaction;
class ReactionSystem;

// Forward declarations of keyword types
namespace kwd {

class EquilibriumPath;
class EquilibriumProblem;
class KineticPath;
class MineralReaction;
class Plot;

} // namespace kwd

/// Initialize a ChemicalPlot object using a kwd::Plot object.
auto initializeChemicalPlot(ChemicalPlot& plot, const kwd::Plot& keyword) -> void;

/// Initialize a MineralReaction object using a kwd::MineralReaction object.
auto initializeMineralReaction(MineralReaction& reaction, const kwd::MineralReaction& keyword) -> void;

/// Initialize an EquilibriumProblem object using a kwd::EquilibriumProblem object.
auto initializeEquilibriumProblem(EquilibriumProblem& problem, const kwd::EquilibriumProblem& keyword) -> void;

/// Initialize an EquilibriumPath object using a kwd::EquilibriumPath object.
auto initializeEquilibriumPath(EquilibriumPath& path, const kwd::EquilibriumPath& keyword) -> void;

/// Initialize a KineticPath object using a kwd::KineticPath object.
auto initializeKineticPath(KineticPath& path, const kwd::KineticPath& keyword) -> void;

} // namespace Reaktoro
