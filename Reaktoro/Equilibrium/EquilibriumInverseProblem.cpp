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

#include "EquilibriumInverseProblem.hpp"

// C++ includes
#include <unordered_map>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ElementUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalProperty.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktoro {
namespace {

/// A type used to define the residual of an equilibrium constraint and its partial derivatives.
struct ResidualEquilibriumConstraint
{
    /// Construct a default EquilibriumConstraintResult instance
    ResidualEquilibriumConstraint()
    {}

    /// Construct an EquilibriumConstraintResult instance from a ChemicalScalar instance
    ResidualEquilibriumConstraint(const ChemicalScalar& scalar)
    : val(scalar.val),
      ddT(scalar.ddT),
      ddP(scalar.ddP),
      ddn(scalar.ddn) {}

    /// The residual value of the equilibrium constraint.
    double val = 0.0;

    /// The partial derivative of the residual w.r.t. temperature.
    double ddT = 0.0;

    /// The partial derivative of the residual w.r.t. pressure.
    double ddP = 0.0;

    /// The partial derivatives of the residual w.r.t. titrant amounts x.
    Vector ddx;

    /// The partial derivatives of the residual w.r.t. species amounts n.
    Vector ddn;
};

/// A type used to represent a titrant and its attributes.
struct Titrant
{
    /// The name of the titrant.
    std::string name;

    /// The elemental formula of the titrant w.r.t. elements in the chemical system.
    Vector formula;

    /// The molar mass of the titrant (in units of kg/mol)
    double molar_mass = 0.0;
};

/// Return the elemental formula of a titrant given as a single compound or as a mixture of compounds.
auto titrantElementalFormula(std::string titrant) -> std::unordered_map<std::string, double>
{
    // First, split on `;` to check if titrant has mixture format like `1 kg H2O; 1 mol NaCl`
    const auto words = split(titrant, ";");

    // Check if size is 1, so that titrant is in a compound-format like `H2O`
    if(words.size() == 1)
        return Reaktoro::elements(split(titrant).back());

    // The elemntal formula of the titrant in solution format
    std::unordered_map<std::string, double> formula;

    // Otherwise, process each word
    for(auto word : words)
    {
        // Split on space ` ` to get three words, like `1 kg H2O` results in ("1", "kg", "H2O")
        const auto triplet = split(word);

        // Define some auxiliary references
        const auto number = tofloat(triplet[0]);
        const auto units = triplet[1];
        const auto name = triplet[2];

        // Get the elements that compose the current titrant component
        const auto elements = Reaktoro::elements(name);

        // Get the molar mass of the current titrant component
        const auto molar_mass = Reaktoro::molarMass(elements);

        // The conversion factor from given units to mol
        double factor = 1;

        // Check if the type of units, either mass or mol
        if(units::convertible(units, "kg"))
            factor = units::convert(1.0, units, "kg")/molar_mass; // convert from mass units to kg, then from kg to mol
        else if(units::convertible(units, "mol"))
            factor = units::convert(1.0, units, "mol"); // convert from molar units to mol
        else RuntimeError("Could not create the titrant with name `" + titrant + "`.",
            "The units of each component in a multi-component titrant must be convertible to mol or kg.");

        // Determine the molar coefficient of the titrant (e.g., 1 kg H2O ~ 55.508 mol H2O)
        const auto coeff = factor * number;

        // Update the elemental formula of the titrant in solution format
        for(auto pair : elements)
            if(formula.count(pair.first))
                formula[pair.first] += coeff * pair.second;
            else
                formula.insert({pair.first, coeff * pair.second});
    }

    return formula;
}

/// Return a titrant created from a given name and a chemical system
auto createTitrant(std::string titrant, const ChemicalSystem& system) -> Titrant
{
    Titrant res;
    res.name = titrant;
    res.formula = zeros(system.numElements());

    Index idx = system.indexSpecies(titrant);

    if(idx < system.numSpecies())
    {
        Species species = system.species(idx);
        res.molar_mass = species.molarMass();
        for(auto pair : species.elements())
            res.formula[system.indexElement(pair.first.name())] = pair.second;
    }
    else
    {
        auto elemental_formula = titrantElementalFormula(titrant);
        for(auto pair : elemental_formula)
            res.formula[system.indexElement(pair.first)] = pair.second;
        res.molar_mass = Reaktoro::molarMass(elemental_formula);
    }

    return res;
}

/// A type used to define the functional signature of an equilibrium constraint.
using EquilibriumConstraint =
    std::function<ResidualEquilibriumConstraint
        (VectorConstRef, const ChemicalState&, const ChemicalProperties&)>;

} // namespace

struct EquilibriumInverseProblem::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The equilibrium problem used to manage temperature, pressure, and mixture of compounds
    EquilibriumProblem problem;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The equilibrium constraint functions
    std::vector<EquilibriumConstraint> constraints;

    /// The unknowns in the problem and their indices
    std::unordered_map<std::string, Index> unknowns;

    /// The names of the titrants
    std::vector<Titrant> titrants;

    /// The formula matrix of the titrants
    Matrix formula_matrix_titrants;

    /// The initial guess of the titrants (in units of mol)
    Vector titrant_initial_amounts;

    /// Construct an Impl instance
    Impl(const ChemicalSystem& system)
    : system(system), problem(system)
    {
        // Initialize a default partition for the chemical system
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        partition = partition_;
        problem.setPartition(partition);
    }

    /// Add a species amount constraint to the inverse equilibrium problem.
    auto addSpeciesAmountConstraint(std::string species, double value) -> void
    {
        // The index of the species
        const Index ispecies = system.indexSpeciesWithError(species);

        // The number of species in the system
        const Index num_species = system.numSpecies();

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar ni(num_species);

        // Set the parial molar derivative of `ni`
        ni.ddn[ispecies] = 1.0;

        // Define the amount constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            ni.val = state.speciesAmount(ispecies);
            return ni - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a species chemical potential constraint to the inverse equilibrium problem.
    auto addSpeciesChemicalPotentialConstraint(std::string species, double value) -> void
    {
        // The index of the species
        const Index ispecies = system.indexSpeciesWithError(species);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar ui;

        // Define the chemical potential constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            ui = properties.chemicalPotentials()[ispecies];
            return ui - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a species activity constraint to the inverse equilibrium problem.
    auto addSpeciesActivityConstraint(std::string species, double value) -> void
    {
        // The index of the species
        const Index ispecies = system.indexSpeciesWithError(species);

        // The ln of the given activity value
        const double ln_val = std::log(value);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar ln_ai;

        // Define the activity constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            ln_ai = properties.lnActivities()[ispecies];
            return ln_ai - ln_val;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a pE constraint to the inverse equilibrium problem.
    auto add_pE_Constraint(double value) -> void
    {
        // The chemical property function that calculates pE
        const auto pE = ChemicalProperty::pE(system);

        // Define the activity constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            return pE(properties) - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a Eh constraint to the inverse equilibrium problem.
    auto addEhConstraint(double value) -> void
    {
        // The chemical property function that calculates Eh
        const auto Eh = ChemicalProperty::Eh(system);

        // Define the activity constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            return Eh(properties) - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a total alkalinity constraint to the inverse equilibrium problem.
    auto addAlkalinityConstraint(double value) -> void
    {
        // The chemical property function that calculates alkalinity
        const auto alk = ChemicalProperty::alkalinity(system);

        // Define the activity constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            return alk(properties) - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    auto addVolumeConstraint(double value) -> void
    {
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            return properties.volume() - value;
        };

        constraints.push_back(f);
    }

    auto addInternalEnergyConstraint(double value) -> void
    {
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            return properties.internalEnergy() - value;
        };

        constraints.push_back(f);
    }

    auto addEnthalpyConstraint(double value) -> void
    {
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            return properties.enthalpy() - value;
        };

        constraints.push_back(f);
    }

    /// Add a phase amount constraint to the inverse equilibrium problem.
    auto addPhaseAmountConstraint(std::string phase, double value) -> void
    {
        // The index of the species
        const Index iphase = system.indexPhaseWithError(phase);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar np;

        // Define the phase volume constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            np = properties.phaseAmounts()[iphase];
            return np - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a phase mass constraint to the inverse equilibrium problem.
    auto addPhaseMassConstraint(std::string phase, double value) -> void
    {
        // The index of the species
        const Index iphase = system.indexPhaseWithError(phase);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar mass;

        // Define the phase volume constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            mass = properties.phaseMasses()[iphase];
            return mass - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a phase volume constraint to the inverse equilibrium problem.
    auto addPhaseVolumeConstraint(std::string phase, double value) -> void
    {
        // The index of the phase
        const Index iphase = system.indexPhaseWithError(phase);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar Vp;

        // Define the phase volume constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            Vp = properties.phaseVolumes()[iphase];
            return Vp - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a sum of phase volumes constraint to the inverse equilibrium problem.
    auto addSumPhaseVolumesConstraint(const std::vector<std::string>& phases, double value) -> void
    {
        // The indices of the phases
        const Indices iphases = system.indicesPhases(phases);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar Vp;

        // Define the phase volume constraint function
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            Vp = sum(rows(properties.phaseVolumes(), iphases));
            return Vp - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Return the index of a titrant.
    auto indexTitrant(std::string titrant) const -> Index
    {
        Index idx = 0;
        for(const auto& t : titrants)
            if(t.name == titrant) return idx; else ++idx;
        return idx;
    }

    /// Set the initial guess of a titrant.
    auto setTitrantInitialAmount(std::string titrant, double amount) -> void
    {
        // Get the index of the titrant
        const Index ititrant = indexTitrant(titrant);

        // Assert this titrant has been added already
        Assert(ititrant < titrants.size(), "Could not set the initial guess "
            "of titrant `" + titrant + "`.", "This titrant has not been added yet.");

        // Set the initial guess of the titrant
        titrant_initial_amounts[ititrant] = amount;
    }

    /// Add an unknown to the inverse equilibrium problem.
    /// @param name The name of the unknown variable (e.g., Temperature, Pressure, HCl, CO2)
    auto addUnknown(std::string name) -> void
    {
        auto i = unknowns.find(name);
        if(i == unknowns.end())
            unknowns.insert({name, unknowns.size()});
    }

    /// Add a titrant to the inverse equilibrium problem.
    auto addTitrant(std::string titrant) -> void
    {
        // Assert that the new titrant has not been added before
        Assert(indexTitrant(titrant) >= titrants.size(),
            "Could not add the titrant " + titrant + " to the inverse problem.",
            "This titrant has been added before.");

        // Update the list of titrants
        titrants.push_back(createTitrant(titrant, system));

        // The number of elements in the system
        const Index num_elements = system.numElements();

        // The index of the last titrant
        const Index ilast = titrants.size() - 1;

        // Update the formula matrix of the formulas
        formula_matrix_titrants.conservativeResize(num_elements, titrants.size());

        // Set the last created column to zero
        formula_matrix_titrants.col(ilast) = titrants.back().formula;

        // Resize the vector of titrant initial guess
        titrant_initial_amounts.conservativeResize(titrants.size());

        // Initialize the initial guess of the new titrant to zero
        titrant_initial_amounts[ilast] = 0.0;
    }

    /// Set two titrants as mutually exclusive.
    auto setAsMutuallyExclusive(std::string titrant1, std::string titrant2) -> void
    {
        // Common error message
        auto errormsg = "Could not set `" + titrant1 + "` and " + titrant2 + "` "
            "as mutually exclusive titrants.";

        // Get the indices of the titrants
        const Index i1 = unknowns.find(titrant1)->second;
        const Index i2 = unknowns.find(titrant2)->second;

        // Check if the two titrants are different
        Assert(titrant1 != titrant2, errormsg,
            "They must have different identifiers.");

        // Check if the first titrant has been added before
        Assert(i1 < unknowns.size(), errormsg,
            "The titrant `" + titrant1 + "` has not been added before.");

        // Check if the second titrant has been added before
        Assert(i2 < unknowns.size(), errormsg,
            "The titrant `" + titrant2 + "` has not been added before.");

        // The smoothing parameter for the mutually exclusive constraint function
        const double tau = 1e-20;

        // Define the mutually exclusive constraint function
        ResidualEquilibriumConstraint res;
        EquilibriumConstraint f = [=](VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) mutable
        {
            const Index Nx = x.rows();

            res.ddx.resize(Nx);

            const double x1 = x[i1];
            const double x2 = x[i2];

            res.val = x1*x2 - tau;
            res.ddx = x1 * unit(Nx, i2) + x2 * unit(Nx, i1);

            return res;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Return the residual of the equilibrium constraints and their partial molar derivatives.
    auto residualEquilibriumConstraints(VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) const -> ResidualEquilibriumConstraints
    {
        const Index num_species = system.numSpecies();
        const Index num_constraints = constraints.size();
        const Index Ee = partition.numEquilibriumElements();
        ResidualEquilibriumConstraints res;
        res.val = zeros(num_constraints);
        res.ddx = zeros(num_constraints, num_constraints);
        res.ddu = zeros(num_constraints, 2 + Ee);
        res.ddn = zeros(num_constraints, num_species);
        ResidualEquilibriumConstraint aux;
        for(Index i = 0; i < num_constraints; ++i)
        {
            aux = constraints[i](x, state, properties);
            res.val[i] = aux.val;
            res.ddu(i, 0) = aux.ddT;
            res.ddu(i, 1) = aux.ddP;
            if(aux.ddx.size())
                res.ddx.row(i) = aux.ddx;
            if(aux.ddn.size())
                res.ddn.row(i) = aux.ddn;
        }
        return res;
    }

    /// Return the initial guess for the unknowns in the inverse chemical equilibrium problem.
    auto initialGuessUnknowns() const -> Vector
    {
        Vector x = zeros(unknowns.size());

        for(auto pair : unknowns)
        {
            if(pair.first == "Temperature")
                x[pair.second] = 0.0; // x[0] = dT = 0 as initial guess
            else if(pair.first == "Pressure")
                x[pair.second] = 0.0; // x[1] = dP = 0 as initial guess
            else {
                const auto i = indexTitrant(pair.first);
                x[pair.second] = titrant_initial_amounts[i];
            }
        }

        return x;
    }

    /// Return the coefficient matrix that relates forward input variables with the unknowns in the problem.
    auto unknownsCoefficientMatrix() const -> Matrix
    {
        const auto iee = partition.indicesEquilibriumElements();
        const auto Ee = partition.numEquilibriumElements();

        Matrix C = zeros(2 + Ee, unknowns.size());

        for(auto pair : unknowns)
        {
            if(pair.first == "Temperature")
                C(0, pair.second) = 1.0;
            else if(pair.first == "Pressure")
                C(1, pair.second) = 1.0;
            else {
                const auto i = indexTitrant(pair.first);
                C.col(pair.second).tail(Ee) = titrants[i].formula(iee);
            }
        }

        return C;
    }
};

EquilibriumInverseProblem::EquilibriumInverseProblem(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumInverseProblem::EquilibriumInverseProblem(const EquilibriumInverseProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumInverseProblem::~EquilibriumInverseProblem()
{}

auto EquilibriumInverseProblem::operator=(EquilibriumInverseProblem other) -> EquilibriumInverseProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumInverseProblem::setPartition(const Partition& partition) -> EquilibriumInverseProblem&
{
    pimpl->setPartition(partition);
    return *this;
}

auto EquilibriumInverseProblem::setTemperature(double val) -> EquilibriumInverseProblem&
{
    pimpl->problem.setTemperature(val);
    return *this;
}

auto EquilibriumInverseProblem::setTemperature(double val, std::string units) -> EquilibriumInverseProblem&
{
    pimpl->problem.setTemperature(val, units);
    return *this;
}

auto EquilibriumInverseProblem::setPressure(double val) -> EquilibriumInverseProblem&
{
    pimpl->problem.setPressure(val);
    return *this;
}

auto EquilibriumInverseProblem::setPressure(double val, std::string units) -> EquilibriumInverseProblem&
{
    pimpl->problem.setPressure(val, units);
    return *this;
}

auto EquilibriumInverseProblem::setElementInitialAmounts(VectorConstRef values) -> EquilibriumInverseProblem&
{
    pimpl->problem.setElementAmounts(values);
    return *this;
}

auto EquilibriumInverseProblem::add(std::string name, double amount, std::string units) -> EquilibriumInverseProblem&
{
    pimpl->problem.add(name, amount, units);
    return *this;
}

auto EquilibriumInverseProblem::fixSpeciesAmount(std::string species, double value, std::string units) -> EquilibriumInverseProblem&
{
    return fixSpeciesAmount(species, value, units, species);
}

auto EquilibriumInverseProblem::fixSpeciesAmount(std::string species, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "mol");
    pimpl->addSpeciesAmountConstraint(species, value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    pimpl->setTitrantInitialAmount(titrant, value);
    return *this;
}

auto EquilibriumInverseProblem::fixSpeciesMass(std::string species, double value, std::string units) -> EquilibriumInverseProblem&
{
    return fixSpeciesMass(species, value, units, species);
}

auto EquilibriumInverseProblem::fixSpeciesMass(std::string species, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    const Index ispecies = pimpl->system.indexSpeciesWithError(species);
    const double molar_mass = pimpl->system.species(ispecies).molarMass();
    value = units::convert(value, units, "kg")/molar_mass;
    return fixSpeciesAmount(species, value, "mol", titrant);
}

auto EquilibriumInverseProblem::fixSpeciesChemicalPotential(std::string species, double value) -> EquilibriumInverseProblem&
{
    pimpl->addSpeciesChemicalPotentialConstraint(species, value);
    return *this;
}

auto EquilibriumInverseProblem::fixSpeciesActivity(std::string species, double value) -> EquilibriumInverseProblem&
{
    pimpl->addSpeciesActivityConstraint(species, value);
    pimpl->addUnknown(species);
    pimpl->addTitrant(species);
    return *this;
}

auto EquilibriumInverseProblem::fixSpeciesActivity(std::string species, double value, std::string titrant) -> EquilibriumInverseProblem&
{
    pimpl->addSpeciesActivityConstraint(species, value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    return *this;
}

auto EquilibriumInverseProblem::fixSpeciesActivity(std::string species, double value, std::string titrant1, std::string titrant2) -> EquilibriumInverseProblem&
{
    pimpl->addSpeciesActivityConstraint(species, value);
    pimpl->addUnknown(titrant1);
    pimpl->addUnknown(titrant2);
    pimpl->addTitrant(titrant1);
    pimpl->addTitrant(titrant2);
    pimpl->setAsMutuallyExclusive(titrant1, titrant2);
    pimpl->setTitrantInitialAmount(titrant1, 1e-6);
    pimpl->setTitrantInitialAmount(titrant2, 1e-6);
    return *this;
}

auto EquilibriumInverseProblem::fixSpeciesFugacity(std::string species, double value, std::string units) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "bar");
    return fixSpeciesActivity(species, value);
}

auto EquilibriumInverseProblem::fixSpeciesFugacity(std::string species, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "bar");
    return fixSpeciesActivity(species, value, titrant);
}

auto EquilibriumInverseProblem::fixPhaseAmount(std::string phase, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "mol");
    pimpl->addPhaseAmountConstraint(phase, value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    pimpl->setTitrantInitialAmount(titrant, value);
    return *this;
}

auto EquilibriumInverseProblem::fixPhaseMass(std::string phase, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "kg");
    pimpl->addPhaseMassConstraint(phase, value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
//    pimpl->setTitrantInitialAmount(titrant, value); // access to molar mass of titrant is needed here for setting adequate initial guess
    return *this;
}

auto EquilibriumInverseProblem::fixPhaseVolume(std::string phase, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "m3");
    pimpl->addPhaseVolumeConstraint(phase, value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    pimpl->setTitrantInitialAmount(titrant, 1e3);
    return *this;
}

auto EquilibriumInverseProblem::fixPhaseSetVolume(const std::vector<std::string>& phases, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "m3");
    pimpl->addSumPhaseVolumesConstraint(phases, value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    pimpl->setTitrantInitialAmount(titrant, 1e3);
    return *this;
}

auto EquilibriumInverseProblem::pH(double value) -> EquilibriumInverseProblem&
{
    const double aHplus = std::pow(10.0, -value);
    return fixSpeciesActivity("H+", aHplus);
}

auto EquilibriumInverseProblem::pH(double value, std::string titrant) -> EquilibriumInverseProblem&
{
    const double aHplus = std::pow(10.0, -value);
    return fixSpeciesActivity("H+", aHplus, titrant);
}

auto EquilibriumInverseProblem::pH(double value, std::string titrant1, std::string titrant2) -> EquilibriumInverseProblem&
{
    const double aHplus = std::pow(10.0, -value);
    return fixSpeciesActivity("H+", aHplus, titrant1, titrant2);
}

auto EquilibriumInverseProblem::pE(double value) -> EquilibriumInverseProblem&
{
    return pE(value, "O2");
}

auto EquilibriumInverseProblem::pE(double value, std::string titrant) -> EquilibriumInverseProblem&
{
    pimpl->add_pE_Constraint(value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    return *this;
}

auto EquilibriumInverseProblem::Eh(double value, std::string units) -> EquilibriumInverseProblem&
{
    return Eh(value, units, "O2");
}

auto EquilibriumInverseProblem::Eh(double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "V");
    pimpl->addEhConstraint(value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    return *this;
}

auto EquilibriumInverseProblem::alkalinity(double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "eq/L");
    pimpl->addAlkalinityConstraint(value);
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
    return *this;
}

auto EquilibriumInverseProblem::fixVolume(double value, std::string units) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "m3");
    pimpl->addVolumeConstraint(value);
    return *this;
}

auto EquilibriumInverseProblem::fixInternalEnergy(double value, std::string units) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "J");
    pimpl->addInternalEnergyConstraint(value);
    return *this;
}

auto EquilibriumInverseProblem::fixEnthalpy(double value, std::string units) -> EquilibriumInverseProblem&
{
    value = units::convert(value, units, "J");
    pimpl->addEnthalpyConstraint(value);
    return *this;
}

auto EquilibriumInverseProblem::unknownTemperature() -> void
{
    pimpl->addUnknown("Temperature");
}

auto EquilibriumInverseProblem::unknownPressure() -> void
{
    pimpl->addUnknown("Pressure");
}

auto EquilibriumInverseProblem::unknownAmountOf(std::string titrant) -> void
{
    pimpl->addUnknown(titrant);
    pimpl->addTitrant(titrant);
}

auto EquilibriumInverseProblem::unknownAmountOfEither(std::string titrant1, std::string titrant2) -> void
{
    pimpl->addUnknown(titrant1);
    pimpl->addUnknown(titrant2);
    pimpl->addTitrant(titrant1);
    pimpl->addTitrant(titrant2);
    pimpl->setAsMutuallyExclusive(titrant1, titrant2);
}

auto EquilibriumInverseProblem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumInverseProblem::partition() const -> const Partition&
{
    return pimpl->partition;
}

auto EquilibriumInverseProblem::temperature() const -> double
{
    return pimpl->problem.temperature();
}

auto EquilibriumInverseProblem::pressure() const -> double
{
    return pimpl->problem.pressure();
}

auto EquilibriumInverseProblem::numConstraints() const -> Index
{
    return pimpl->constraints.size();
}

auto EquilibriumInverseProblem::numUnknowns() const -> Index
{
    return pimpl->unknowns.size();
}

auto EquilibriumInverseProblem::numTitrants() const -> Index
{
    return pimpl->titrants.size();
}

auto EquilibriumInverseProblem::formulaMatrixTitrants() const -> Matrix
{
    return pimpl->formula_matrix_titrants;
}

auto EquilibriumInverseProblem::unknownsCoefficientMatrix() const -> Matrix
{
    return pimpl->unknownsCoefficientMatrix();
}

auto EquilibriumInverseProblem::elementInitialAmounts() const -> Vector
{
    return pimpl->problem.elementAmounts();
}

auto EquilibriumInverseProblem::titrantInitialAmounts() const -> Vector
{
    return pimpl->titrant_initial_amounts;
}

auto EquilibriumInverseProblem::residualEquilibriumConstraints(VectorConstRef x, const ChemicalState& state, const ChemicalProperties& properties) const -> ResidualEquilibriumConstraints
{
    return pimpl->residualEquilibriumConstraints(x, state, properties);
}

auto EquilibriumInverseProblem::initialGuessUnknowns() const -> Vector
{
    return pimpl->initialGuessUnknowns();
}

} // namespace Reaktoro
