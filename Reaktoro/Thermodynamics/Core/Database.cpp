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

#include "Database.hpp"

// C++ includes
#include <clocale>
#include <map>
#include <set>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/GlobalOptions.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Databases/DatabaseUtils.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/FluidSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

// miniz includes
#include <miniz/zip_file.hpp>

// pugixml includes
#include <pugixml.hpp>
using namespace pugi;

namespace Reaktoro {
namespace {

/// Auxiliary types for a map of aqueous, fluid, and mineral species
using ElementMap        = std::map<std::string, Element>;
using AqueousSpeciesMap = std::map<std::string, AqueousSpecies>;
using FluidSpeciesMap   = std::map<std::string, FluidSpecies>;
using MineralSpeciesMap = std::map<std::string, MineralSpecies>;

auto errorNonExistentSpecies(std::string type, std::string name) -> void
{
    Exception exception;
    exception.error << "Cannot get an instance of the " << type << " species `" << name << "` in the database.";
    exception.reason << "There is no such species in the database.";
    RaiseError(exception);
}

auto parseDissociation(std::string dissociation) -> std::map<std::string, double>
{
    std::map<std::string, double> equation;
    auto words = split(dissociation, " ");
    for(const auto& word : words)
    {
        auto pair = split(word, ":");
        equation.emplace(pair[1], tofloat(pair[0]));
    }
    return equation;
}

auto parseReactionInterpolatedThermoProperties(const xml_node& node) -> ReactionThermoInterpolatedProperties
{
    // Get the data values of the children nodes
    std::vector<double> temperatures     = tofloats(node.child("Temperatures").text().get());
    std::vector<double> pressures        = tofloats(node.child("Pressures").text().get());
    std::vector<double> pk               = tofloats(node.child("pk").text().get());
    std::vector<double> lnk              = tofloats(node.child("lnk").text().get());
    std::vector<double> logk             = tofloats(node.child("logk").text().get());
    std::vector<double> gibbs_energy     = tofloats(node.child("G").text().get());
    std::vector<double> helmholtz_energy = tofloats(node.child("A").text().get());
    std::vector<double> internal_energy  = tofloats(node.child("U").text().get());
    std::vector<double> enthalpy         = tofloats(node.child("H").text().get());
    std::vector<double> entropy          = tofloats(node.child("S").text().get());
    std::vector<double> volume           = tofloats(node.child("V").text().get());
    std::vector<double> heat_capacity_cp = tofloats(node.child("Cp").text().get());
    std::vector<double> heat_capacity_cv = tofloats(node.child("Cv").text().get());

    // Convert `pk` to `lnk`, where `pk = -log(k) = -ln(k)/ln(10)`
    const double ln_10 = std::log(10.0);
    if(!pk.empty() && lnk.empty())
    {
        lnk.resize(pk.size());
        for(unsigned i = 0; i < pk.size(); ++i)
            lnk[i] = -pk[i] * ln_10;
    }

    // Convert `logk` to `lnk`, where `log(k) = ln(k)/ln(10)`
    if(!logk.empty() && lnk.empty())
    {
        lnk.resize(logk.size());
        for(unsigned i = 0; i < logk.size(); ++i)
            lnk[i] = logk[i] * ln_10;
    }

    // Get the temperature and pressure units
    std::string tunits = node.child("Temperatures").attribute("units").as_string();
    std::string punits = node.child("Pressures").attribute("units").as_string();

    // Get the names and stoichiometries of the species that define the reaction
    std::string equation = node.child("Equation").text().get();

    // Check if element `temperatures` was provided, if not set default to 25 celsius
    if(temperatures.empty()) temperatures.push_back(25.0);

    // Check if element `pressures` was provided, if not set default to 1 bar
    if(pressures.empty()) pressures.push_back(1.0);

    // Check if temperature units was provided, if not set default to celsius
    if(tunits.empty()) tunits = "celsius";

    // Check if pressure units was provided, if not set default to bar
    if(punits.empty()) punits = "bar";

    // Convert temperatures and pressures to standard units (kelvin and pascal respectively)
    for(auto& x : temperatures) x = units::convert(x, tunits, "kelvin");
    for(auto& x : pressures)    x = units::convert(x, punits, "pascal");

    // Define a lambda function to generate a bilinear interpolator from a vector
    auto bilinear_interpolator = [&](const std::vector<double>& data) -> BilinearInterpolator
    {
        if(data.empty())  return BilinearInterpolator();
        return BilinearInterpolator(temperatures, pressures, data);
    };

    // Define a lambda function to generate a bilinear interpolator of the Gibbs energy of a reaction from its lnk
    auto gibbs_energy_from_lnk = [](const BilinearInterpolator& lnk)
    {
        const double R = universalGasConstant;
        auto f = [=](double T, double P) { return -R*T*lnk(T, P); };
        return BilinearInterpolator(lnk.xCoordinates(), lnk.yCoordinates(), f);
    };

    // Initialize the properties thermodynamic properties of the reaction
    ReactionThermoInterpolatedProperties data;
    data.equation         = equation;
    data.lnk              = bilinear_interpolator(lnk);
    data.gibbs_energy     = gibbs_energy.empty() ? gibbs_energy_from_lnk(data.lnk) : bilinear_interpolator(gibbs_energy);
    data.helmholtz_energy = bilinear_interpolator(helmholtz_energy);
    data.internal_energy  = bilinear_interpolator(internal_energy);
    data.enthalpy         = bilinear_interpolator(enthalpy);
    data.entropy          = bilinear_interpolator(entropy);
    data.volume           = bilinear_interpolator(volume);
    data.heat_capacity_cp = bilinear_interpolator(heat_capacity_cp);
    data.heat_capacity_cv = bilinear_interpolator(heat_capacity_cv);

    return data;
}

auto parseSpeciesInterpolatedThermoProperties(const xml_node& node) -> SpeciesThermoInterpolatedProperties
{
    // Get the data values of the children nodes
    std::vector<double> temperatures     = tofloats(node.child("Temperatures").text().get());
    std::vector<double> pressures        = tofloats(node.child("Pressures").text().get());
    std::vector<double> gibbs_energy     = tofloats(node.child("G").text().get());
    std::vector<double> helmholtz_energy = tofloats(node.child("A").text().get());
    std::vector<double> internal_energy  = tofloats(node.child("U").text().get());
    std::vector<double> enthalpy         = tofloats(node.child("H").text().get());
    std::vector<double> entropy          = tofloats(node.child("S").text().get());
    std::vector<double> volume           = tofloats(node.child("V").text().get());
    std::vector<double> heat_capacity_cp = tofloats(node.child("Cp").text().get());
    std::vector<double> heat_capacity_cv = tofloats(node.child("Cv").text().get());

    // Get the temperature and pressure units
    std::string tunits = node.child("Temperatures").attribute("units").as_string();
    std::string punits = node.child("Pressures").attribute("units").as_string();

    // Check if element `temperatures` was provided, if not set default to 25 celsius
    if(temperatures.empty()) temperatures.push_back(25.0);

    // Check if element `pressures` was provided, if not set default to 1 bar
    if(pressures.empty()) pressures.push_back(1.0);

    // Check if temperature units was provided, if not set default to celsius
    if(tunits.empty()) tunits = "celsius";

    // Check if pressure units was provided, if not set default to bar
    if(punits.empty()) punits = "bar";

    // Convert temperatures and pressures to standard units
    for(auto& x : temperatures) x = units::convert(x, tunits, "kelvin");
    for(auto& x : pressures)    x = units::convert(x, punits, "pascal");

    // Define a lambda function to generate a bilinear interpolator from a vector
    auto bilinear_interpolator = [&](const std::vector<double>& data) -> BilinearInterpolator
    {
        if(data.empty())  return BilinearInterpolator();
        return BilinearInterpolator(temperatures, pressures, data);
    };

    // Initialize the properties thermodynamic properties of the species
    SpeciesThermoInterpolatedProperties data;
    data.gibbs_energy     = bilinear_interpolator(gibbs_energy);
    data.helmholtz_energy = bilinear_interpolator(helmholtz_energy);
    data.internal_energy  = bilinear_interpolator(internal_energy);
    data.enthalpy         = bilinear_interpolator(enthalpy);
    data.entropy          = bilinear_interpolator(entropy);
    data.volume           = bilinear_interpolator(volume);
    data.heat_capacity_cp = bilinear_interpolator(heat_capacity_cp);
    data.heat_capacity_cv = bilinear_interpolator(heat_capacity_cv);

    return data;
}

auto as_int(const xml_node& node, const char* childname, int if_empty=std::numeric_limits<int>::infinity()) -> int
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().as_int();
}

auto as_double(const xml_node& node, const char* childname, double if_empty=std::numeric_limits<double>::infinity()) -> double
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().as_double();
}

auto parseAqueousSpeciesThermoParamsHKF(const xml_node& node) -> std::optional<AqueousSpeciesThermoParamsHKF>
{
    AqueousSpeciesThermoParamsHKF hkf;
    hkf.Gf   = as_double(node, "Gf");
    hkf.Hf   = as_double(node, "Hf");
    hkf.Sr   = as_double(node, "Sr");
    hkf.a1   = as_double(node, "a1");
    hkf.a2   = as_double(node, "a2");
    hkf.a3   = as_double(node, "a3");
    hkf.a4   = as_double(node, "a4");
    hkf.c1   = as_double(node, "c1");
    hkf.c2   = as_double(node, "c2");
    hkf.wref = as_double(node, "wref");
    return hkf;
}

auto parseFluidSpeciesThermoParamsHKF(const xml_node& node)->std::optional<FluidSpeciesThermoParamsHKF>
{
    FluidSpeciesThermoParamsHKF hkf;
    hkf.Gf   = as_double(node, "Gf");
    hkf.Hf   = as_double(node, "Hf");
    hkf.Sr   = as_double(node, "Sr");
    hkf.a    = as_double(node, "a");
    hkf.b    = as_double(node, "b");
    hkf.c    = as_double(node, "c");
    hkf.Tmax = as_double(node, "Tmax");
    return hkf;
}

auto parseMineralSpeciesThermoParamsHKF(const xml_node& node) -> std::optional<MineralSpeciesThermoParamsHKF>
{
    MineralSpeciesThermoParamsHKF hkf;
    hkf.Gf      = as_double(node, "Gf");
    hkf.Hf      = as_double(node, "Hf");
    hkf.Sr      = as_double(node, "Sr");
    hkf.Vr      = as_double(node, "Vr");
    hkf.nptrans = as_int(node, "NumPhaseTrans");
    hkf.Tmax    = as_double(node, "Tmax");

    if(hkf.nptrans == 0)
    {
        hkf.a.push_back(as_double(node, "a"));
        hkf.b.push_back(as_double(node, "b"));
        hkf.c.push_back(as_double(node, "c"));
    }
    else
    {
        for(int i = 0; i <= hkf.nptrans; ++i)
        {
            std::stringstream str;
            str << "TemperatureRange" << i;

            auto temperature_range = node.child(str.str().c_str());

            hkf.a.push_back(as_double(temperature_range, "a"));
            hkf.b.push_back(as_double(temperature_range, "b"));
            hkf.c.push_back(as_double(temperature_range, "c"));

            if(i < hkf.nptrans)
            {
                hkf.Ttr.push_back(as_double(temperature_range, "Ttr"));

                // Set zero the non-available transition values
                hkf.Htr.push_back(as_double(temperature_range, "Htr", 0.0));
                hkf.Vtr.push_back(as_double(temperature_range, "Vtr", 0.0));
                hkf.dPdTtr.push_back(as_double(temperature_range, "dPdTtr", 0.0));
            }
        }
    }

    return hkf;
}

auto parseAqueousSpeciesThermoData(const xml_node& node) -> AqueousSpeciesThermoData
{
    AqueousSpeciesThermoData thermo;

    if(!node.child("Properties").empty())
        thermo.properties = parseSpeciesInterpolatedThermoProperties(node.child("Properties"));

    if(!node.child("Reaction").empty())
        thermo.reaction = parseReactionInterpolatedThermoProperties(node.child("Reaction"));

    if(!node.child("HKF").empty())
        thermo.hkf = parseAqueousSpeciesThermoParamsHKF(node.child("HKF"));

    return thermo;
}

auto parseFluidSpeciesThermoData(const xml_node& node) -> FluidSpeciesThermoData
{
    FluidSpeciesThermoData thermo;

    if(!node.child("Properties").empty())
        thermo.properties = parseSpeciesInterpolatedThermoProperties(node.child("Properties"));

    if(!node.child("Reaction").empty())
        thermo.reaction = parseReactionInterpolatedThermoProperties(node.child("Reaction"));

    if(!node.child("HKF").empty())
        thermo.hkf = parseFluidSpeciesThermoParamsHKF(node.child("HKF"));

    return thermo;
}

auto parseMineralSpeciesThermoData(const xml_node& node) -> MineralSpeciesThermoData
{
    MineralSpeciesThermoData thermo;

    if(!node.child("Properties").empty())
        thermo.properties = parseSpeciesInterpolatedThermoProperties(node.child("Properties"));

    if(!node.child("Reaction").empty())
        thermo.reaction = parseReactionInterpolatedThermoProperties(node.child("Reaction"));

    if(!node.child("HKF").empty())
        thermo.hkf = parseMineralSpeciesThermoParamsHKF(node.child("HKF"));

    return thermo;
}

template<typename SpeciesType, typename SpeciesFunction>
auto collectSpecies(const std::map<std::string, SpeciesType>& map, const SpeciesFunction& fn) -> std::vector<SpeciesType>
{
    std::set<SpeciesType> species;
    for(const auto& entry : map)
        if(fn(entry.second))
            species.insert(entry.second);
    return std::vector<SpeciesType>(species.begin(), species.end());
}

template<typename SpeciesType>
auto speciesWithElements(const std::vector<std::string>& elements, const std::map<std::string, SpeciesType>& map) -> std::vector<SpeciesType>
{
    auto f = [&](const SpeciesType& species)
    {
        // Check if the given species has all chemical elements in the list of elements (ignore charge Z)
        for(auto pair : species.elements())
            if(pair.first.name() != "Z" && !contained(pair.first.name(), elements))
                return false;
        return true;
    };

    return collectSpecies(map, f);
}

} // namespace

//A guard object to guarantee the return of original locale
class ChangeLocale
{
public:
    ChangeLocale() = delete;
    explicit ChangeLocale(const char* new_locale) : old_locale(std::setlocale(LC_NUMERIC, nullptr))
    {
        std::setlocale(LC_NUMERIC, new_locale);
    }
    ~ChangeLocale()
    {
        std::setlocale(LC_NUMERIC, old_locale.c_str());
    }

private:
    const std::string old_locale;

};

struct Database::Impl
{
    /// The set of all elements in the database
    ElementMap element_map;

    /// The set of all aqueous species in the database
    AqueousSpeciesMap aqueous_species_map;

    /// The set of all gaseous species in the database
    FluidSpeciesMap gaseous_species_map;

    /// The set of all liquid species in the database
    FluidSpeciesMap liquid_species_map;

    /// The set of all fluid species in the database
    FluidSpeciesMap fluid_species_map;

    /// The set of all mineral species in the database
    MineralSpeciesMap mineral_species_map;

    Impl() = default;

    Impl(std::string filename)
    {
        const auto guard = ChangeLocale("C");

        // Create the XML document
        xml_document doc;

        // Load the xml database file
        auto result = doc.load_file(filename.c_str());

        // Check if result is not ok, and then try a built-in database with same name
        if(!result)
        {
            // Search for a built-in database
            std::string builtin = database(filename);

            // If not empty, use the built-in database to create the xml doc
            if(!builtin.empty()) result = doc.load_string(builtin.c_str());
        }

        // Ensure either a database file path was correctly given, or a built-in database
        if(!result)
        {
            std::string names;
            for(auto const& s : databases()) names += s + " ";
            RuntimeError("Could not initialize the Database instance with given database name `" + filename + "`.",
                "This name either points to a non-existent database file, or it is not one of the "
                "built-in database files in Reaktoro. The built-in databases are: " + names + ".");
        }

        // Parse the xml document
        parse(doc, filename);
    }

    template<typename Key, typename Value>
    auto collectValues(const std::map<Key, Value>& map) -> std::vector<Value>
    {
        std::vector<Value> species;
        species.reserve(map.size());
        for(const auto& pair : map)
            species.push_back(pair.second);
        return species;
    }

    auto addElement(const Element& element) -> void
    {
        element_map.insert({element.name(), element});
    }

    auto addAqueousSpecies(const AqueousSpecies& species) -> void
    {
        aqueous_species_map.insert({species.name(), species});
    }

    auto addFluidSpecies(const FluidSpecies& species) -> void
    {
        fluid_species_map.insert({ species.name(), species });
    }

    auto addGaseousSpecies(const FluidSpecies& species) -> void
    {
        gaseous_species_map.insert({ species.name(), species });
    }

    auto addLiquidSpecies(const FluidSpecies& species) -> void
    {
        liquid_species_map.insert({ species.name(), species });
    }

    auto addMineralSpecies(const MineralSpecies& species) -> void
    {
        mineral_species_map.insert({species.name(), species});
    }

    auto elements() const-> std::vector<Element>
    {
        std::vector<Element> elements{};
        elements.reserve(element_map.size());
        for(const auto& element : element_map) {
            auto element_copy = Element();
            element_copy.setName(element.second.name());
            element_copy.setMolarMass(element.second.molarMass());
            elements.push_back(element_copy);
        }

        return elements;
    }

    auto aqueousSpecies() -> std::vector<AqueousSpecies>
    {
        return collectValues(aqueous_species_map);
    }

    auto aqueousSpecies(std::string name) const -> const AqueousSpecies&
    {
        if(aqueous_species_map.count(name) == 0)
            errorNonExistentSpecies("aqueous", name);

        return aqueous_species_map.find(name)->second;
    }

    auto fluidSpecies() -> std::vector<FluidSpecies>
    {
        return collectValues(fluid_species_map);
    }

    auto fluidSpecies(std::string name) const -> const FluidSpecies&
    {
        if (fluid_species_map.count(name) == 0)
            errorNonExistentSpecies("fluid", name);

        return fluid_species_map.find(name)->second;
    }

    auto gaseousSpecies() -> std::vector<FluidSpecies>
    {
        return collectValues(gaseous_species_map);
    }

    auto gaseousSpecies(std::string name) const -> const FluidSpecies&
    {
        if(gaseous_species_map.count(name) == 0)
            errorNonExistentSpecies("gaseous", name);

        return gaseous_species_map.find(name)->second;
    }

    auto liquidSpecies() -> std::vector<FluidSpecies>
    {
        return collectValues(liquid_species_map);
    }

    auto liquidSpecies(std::string name) const -> const FluidSpecies&
    {
        if(liquid_species_map.count(name) == 0)
            errorNonExistentSpecies("liquid", name);

        return liquid_species_map.find(name)->second;
    }

    auto mineralSpecies() -> std::vector<MineralSpecies>
    {
        return collectValues(mineral_species_map);
    }

    auto mineralSpecies(std::string name) const -> const MineralSpecies&
    {
        if(mineral_species_map.count(name) == 0)
            errorNonExistentSpecies("mineral", name);

        return mineral_species_map.find(name)->second;
    }

    auto containsAqueousSpecies(std::string species) const -> bool
    {
        return aqueous_species_map.count(species) != 0;
    }

    auto containsFluidSpecies(std::string species) const -> bool
    {
        return fluid_species_map.count(species) != 0;
    }

    auto containsGaseousSpecies(std::string species) const -> bool
    {
        return gaseous_species_map.count(species) != 0;
    }

    auto containsLiquidSpecies(std::string species) const -> bool
    {
        return liquid_species_map.count(species) != 0;
    }

    auto containsMineralSpecies(std::string species) const -> bool
    {
        return mineral_species_map.count(species) != 0;
    }

    auto aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<AqueousSpecies>
    {
        return speciesWithElements(elements, aqueous_species_map);
    }

    auto fluidSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<FluidSpecies>
    {
        return speciesWithElements(elements, fluid_species_map);
    }

    auto gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<FluidSpecies>
    {
        return speciesWithElements(elements, gaseous_species_map);
    }

    auto liquidSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<FluidSpecies>
    {
        return speciesWithElements(elements, liquid_species_map);
    }

    auto mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<MineralSpecies>
    {
        return speciesWithElements(elements, mineral_species_map);
    }

    auto parse(const xml_document& doc, std::string databasename) -> void
    {
        // Access the database node of the database file
        xml_node database = doc.child("Database");

        // Read all elements in the database
        for(xml_node node : database.children("Element"))
        {
            Element element = parseElement(node);
            element_map[element.name()] = element;
        }

        // Add charge element
        element_map["Z"] = Element();
        element_map["Z"].setName("Z");

        // Read all species in the database
        for(xml_node node : database.children("Species"))
        {
            std::string type = node.child("Type").text().get();
            std::string name = node.child("Name").text().get();
            if (type == "Gaseous")
            {
                FluidSpecies fluid_species = parseFluidSpecies(node);
                fluid_species.setName(name.substr(0, name.size() - 3));
                FluidSpecies gaseous_species = parseFluidSpecies(node);
                gaseous_species.setName(name);
                FluidSpecies liquid_species = parseFluidSpecies(node);
                liquid_species.setName(name.substr(0, name.size() - 3) + "(liq)");
                if (valid(fluid_species))
                {
                    fluid_species_map[fluid_species.name()] = fluid_species;
                    gaseous_species_map[gaseous_species.name()] = gaseous_species;
                    liquid_species_map[liquid_species.name()] = liquid_species;
                }
            }
            else if (type == "Aqueous")
            {
                AqueousSpecies species = parseAqueousSpecies(node);
                if(valid(species))
                    aqueous_species_map[species.name()] = species;
            }
            else if(type == "Mineral")
            {
                MineralSpecies species = parseMineralSpecies(node);
                if(valid(species))
                    mineral_species_map[species.name()] = species;
            }
            else RuntimeError("Could not parse the species `" +
                name + "` with type `" + type + "` in the database `" +
                databasename + "`.", "The type of the species is unknown.");
        }
    }

    auto parseElement(const xml_node& node) -> Element
    {
        Element element;
        element.setName(node.child("Name").text().get());
        element.setMolarMass(node.child("MolarMass").text().as_double() * 1e-3); // convert from g/mol to kg/mol
        return element;
    }

    auto parseElementalFormula(const xml_node& node) -> std::map<Element, double>
    {
        std::string formula = node.child("Elements").text().get();
        std::map<Element, double> elements;
        auto words = split(formula, "()");
        for(unsigned i = 0; i < words.size(); i += 2)
        {
            Assert(element_map.count(words[i]),
                "Cannot parse the elemental formula `" + formula + "`.",
                "The element `" + words[i] + "` is not in the database.");
            elements.emplace(element_map.at(words[i]), tofloat(words[i + 1]));
        }
        if(!node.child("Charge").empty())
        {
            double charge = node.child("Charge").text().as_double();
            elements.emplace(element_map.at("Z"), charge);
        }
        return elements;
    }

    auto parseSpecies(const xml_node& node) -> Species
    {
        // The species instance
        Species species;

        // Set the name of the species
        species.setName(node.child("Name").text().get());

        // Set the chemical formula of the species
        species.setFormula(node.child("Formula").text().get());

        // Set the elements of the species
        species.setElements(parseElementalFormula(node));

        return species;
    }

    auto parseAqueousSpecies(const xml_node& node) -> AqueousSpecies
    {
        // The aqueous species instance
        AqueousSpecies species = parseSpecies(node);

        // Set the elemental charge of the species
        species.setCharge(node.child("Charge").text().as_double());

        // Parse the complex formula of the aqueous species (if any)
        species.setDissociation(parseDissociation(node.child("Dissociation").text().get()));

        // Parse the thermodynamic data of the aqueous species
        species.setThermoData(parseAqueousSpeciesThermoData(node.child("Thermo")));

        return species;
    }
    
    auto parseFluidSpecies(const xml_node& node) -> FluidSpecies
    {
        // The gaseous species instance
        FluidSpecies species = parseSpecies(node);

        // Set the critical temperature of the fluid (gaseous or liquid) species (in units of K)
        if(!node.child("CriticalTemperature").empty())
            species.setCriticalTemperature(node.child("CriticalTemperature").text().as_double());

        // Set the critical pressure of the fluid (gaseous or liquid) species (in units of Pa)
        if(!node.child("CriticalPressure").empty())
            species.setCriticalPressure(node.child("CriticalPressure").text().as_double() * 1e5); // convert from bar to Pa

        // Set the acentric factor of the fluid (gaseous or liquid) species
        if(!node.child("AcentricFactor").empty())
            species.setAcentricFactor(node.child("AcentricFactor").text().as_double());

        // Parse the thermodynamic data of the fluid (gaseous or liquid) species
        species.setThermoData(parseFluidSpeciesThermoData(node.child("Thermo")));

        return species;
    }

    auto parseMineralSpecies(const xml_node& node) -> MineralSpecies
    {
        // The mineral species instance
        MineralSpecies species = parseSpecies(node);

        // Parse the thermodynamic data of the mineral species
        species.setThermoData(parseMineralSpeciesThermoData(node.child("Thermo")));

        return species;
    }

    /// Return true if a species instance has correct and complete data
    template<typename SpeciesType>
    auto valid(const SpeciesType& species) const -> bool
    {
        // Skip validation if even species with missing data should be considered
        if(!global::options.database.exclude_species_with_missing_data)
            return true;

        // Check if species data is available
        if(species.name().empty())
            return false;
        if(species.elements().empty())
            return false;
        if(species.molarMass() <= 0.0 || !std::isfinite(species.molarMass()))
            return false;
        if(species.formula().empty())
            return false;

        // Check if HKF parameters exist, but they are incomplete
        const auto& hkf = species.thermoData().hkf;
        if(hkf && !std::isfinite(hkf.value().Gf))
            return false;
        if(hkf && !std::isfinite(hkf.value().Hf))
            return false;

        return true;
    }
};

Database::Database()
: pimpl(new Impl())
{}

Database::Database(std::string filename)
: pimpl(new Impl(filename))
{}

auto Database::addElement(const Element& element) -> void
{
    pimpl->addElement(element);
}

auto Database::addAqueousSpecies(const AqueousSpecies& species) -> void
{
    pimpl->addAqueousSpecies(species);
}

auto Database::addGaseousSpecies(const FluidSpecies& species) -> void
{
    pimpl->addGaseousSpecies(species);
}

auto Database::addLiquidSpecies(const FluidSpecies& species) -> void
{
    pimpl->addLiquidSpecies(species);
}

auto Database::addFluidSpecies(const FluidSpecies& species) -> void
{
    pimpl->addFluidSpecies(species);
}

auto Database::addMineralSpecies(const MineralSpecies& species) -> void
{
    pimpl->addMineralSpecies(species);
}

auto Database::elements() const -> std::vector<Element>
{
    return pimpl->elements();
}

auto Database::aqueousSpecies() -> std::vector<AqueousSpecies>
{
    return pimpl->aqueousSpecies();
}

auto Database::aqueousSpecies(std::string name) const -> const AqueousSpecies&
{
    return pimpl->aqueousSpecies(name);
}

auto Database::fluidSpecies() -> std::vector<FluidSpecies>
{
    return pimpl->fluidSpecies();
}

auto Database::fluidSpecies(std::string name) const -> const FluidSpecies&
{
    return pimpl->fluidSpecies(name);
}

auto Database::gaseousSpecies() -> std::vector<FluidSpecies>
{
    return pimpl->gaseousSpecies();
}

auto Database::gaseousSpecies(std::string name) const -> const FluidSpecies&
{
    return pimpl->gaseousSpecies(name);
}

auto Database::liquidSpecies() -> std::vector<FluidSpecies>
{
    return pimpl->liquidSpecies();
}

auto Database::liquidSpecies(std::string name) const -> const FluidSpecies&
{
    return pimpl->liquidSpecies(name);
}

auto Database::mineralSpecies() -> std::vector<MineralSpecies>
{
    return pimpl->mineralSpecies();
}

auto Database::mineralSpecies(std::string name) const -> const MineralSpecies&
{
    return pimpl->mineralSpecies(name);
}

auto Database::containsAqueousSpecies(std::string species) const -> bool
{
    return pimpl->containsAqueousSpecies(species);
}

auto Database::containsGaseousSpecies(std::string species) const -> bool
{
    return pimpl->containsGaseousSpecies(species);
}

auto Database::containsLiquidSpecies(std::string species) const -> bool
{
    return pimpl->containsLiquidSpecies(species);
}

auto Database::containsFluidSpecies(std::string species) const -> bool
{
    return pimpl->containsFluidSpecies(species);
}

auto Database::containsMineralSpecies(std::string species) const -> bool
{
    return pimpl->containsMineralSpecies(species);
}

auto Database::aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<AqueousSpecies>
{
    return pimpl->aqueousSpeciesWithElements(elements);
}

auto Database::gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<FluidSpecies>
{
    return pimpl->gaseousSpeciesWithElements(elements);
}

auto Database::liquidSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<FluidSpecies>
{
    return pimpl->liquidSpeciesWithElements(elements);
}

auto Database::fluidSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<FluidSpecies>
{
    return pimpl->fluidSpeciesWithElements(elements);
}

auto Database::mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<MineralSpecies>
{
    return pimpl->mineralSpeciesWithElements(elements);
}

} // namespace Reaktoro

