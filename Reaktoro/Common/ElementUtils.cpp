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

#include "ElementUtils.hpp"

// C++ includes
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>
using std::string;
using std::pair;
using std::unordered_map;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace internal {

/// The list of all known chemical elements
const std::vector<std::string> elements =
{
    "H",  "He", "Li", "Be", "B",   "C",   "N",   "O",   "F",   "Ne",  "Na",  "Mg", "Al", "Si", "P",
    "S",  "Cl", "Ar", "K",  "Ca",  "Sc",  "Ti",  "V",   "Cr",  "Mn",  "Fe",  "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br",  "Kr",  "Rb",  "Sr",  "Y",   "Zr",  "Nb",  "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn",  "Sb",  "Te",  "I",   "Xe",  "Cs",  "Ba",  "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb",  "Dy",  "Ho",  "Er",  "Tm",  "Yb",  "Lu",  "Hf", "Ta", "W",  "Re",
    "Os", "Ir", "Pt", "Au", "Hg",  "Tl",  "Pb",  "Bi",  "Po",  "At",  "Rn",  "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am",  "Cm",  "Bk",  "Cf",  "Es",  "Fm",  "Md",  "No", "Lr", "Rf", "Db",
    "Sg", "Hs", "Bh", "Mt", "Uun", "Uuu", "Uub", "Uut", "Uuq", "Uup", "Uuh"
};

/// The atomic weights of all known chemical elements (in units of g/mol)
const std::unordered_map<std::string, double> atomicWeights =
{
      {"H",   1.00794},    {"He",  4.002602}, {"Li",  6.941},     {"Be",  9.012182},  {"B",  10.811},    {"C",   12.0107}, {"N",   14.0067},   {"O",   15.9994},
      {"F",   18.9984032}, {"Ne",  20.1797},  {"Na",  22.98977},  {"Mg",  24.305},    {"Al", 26.981538}, {"Si",  28.0855}, {"P",   30.973761}, {"S",   32.065},
      {"Cl",  35.453},     {"Ar",  39.948},   {"K",   39.0983},   {"Ca",  40.078},    {"Sc", 44.95591},  {"Ti",  47.867},  {"V",   50.9415},   {"Cr",  51.9961},
      {"Mn",  54.938049},  {"Fe",  55.845},   {"Co",  58.9332},   {"Ni",  58.6934},   {"Cu", 63.546},    {"Zn",  65.409},  {"Ga",  69.723},    {"Ge",  72.64},
      {"As",  74.9216},    {"Se",  78.96},    {"Br",  79.904},    {"Kr",  83.798},    {"Rb", 85.4678},   {"Sr",  87.62},   {"Y",   88.90585},  {"Zr",  91.224},
      {"Nb",  92.90638},   {"Mo",  95.94},    {"Tc",  98.0},      {"Ru",  101.07},    {"Rh", 102.9055},  {"Pd",  106.42},  {"Ag",  107.8682},  {"Cd",  112.411},
      {"In",  114.818},    {"Sn",  118.71},   {"Sb",  121.76},    {"Te",  127.6},     {"I",  126.90447}, {"Xe",  131.293}, {"Cs",  132.90545}, {"Ba",  137.327},
      {"La",  138.9055},   {"Ce",  140.116},  {"Pr",  140.90765}, {"Nd",  144.24},    {"Pm", 145.0},     {"Sm",  150.36},  {"Eu",  151.964},   {"Gd",  157.25},
      {"Tb",  158.92534},  {"Dy",  162.5},    {"Ho",  164.93032}, {"Er",  167.259},   {"Tm", 168.93421}, {"Yb",  173.04},  {"Lu",  174.967},   {"Hf",  178.49},
      {"Ta",  180.9479},   {"W",   183.84},   {"Re",  186.207},   {"Os",  190.23},    {"Ir", 192.217},   {"Pt",  195.078}, {"Au",  196.96655}, {"Hg",  200.59},
      {"Tl",  204.3833},   {"Pb",  207.2},    {"Bi",  208.98038}, {"Po",  209.0},     {"At", 210.0},     {"Rn",  222.0},   {"Fr",  223.0},     {"Ra",  226.0},
      {"Ac",  227.0},      {"Th",  232.0381}, {"Pa",  231.03588}, {"U",   238.02891}, {"Np", 237.0},     {"Pu",  244.0},   {"Am",  243.0},     {"Cm",  247.0},
      {"Bk",  247.0},      {"Cf",  251.0},    {"Es",  252.0},     {"Fm",  257.0},     {"Md", 258.0},     {"No",  259.0},   {"Lr",  262.0},     {"Rf",  261.0},
      {"Db",  262.0},      {"Sg",  266.0},    {"Hs",  264.0},     {"Bh",  277.0},     {"Mt", 268.0},     {"Uun", 281.0},   {"Uuu", 272.0},     {"Uub", 285.0},
      {"Uut", 284.0},      {"Uuq", 289.0},    {"Uup", 288.0},     {"Uuh", 292.0}
};

auto parseElementAtom(string::iterator begin, string::iterator end) -> pair<string, string::iterator>
{
    assert(isupper(*begin) && "***Error*** `parseElementAtom` requires the first character to be uppercase.");
    if(begin == end) return {"", begin};
    auto endelement = std::find_if(begin + 1, end, [](char c){return isupper(c) || !isalpha(c);});
    string element = string(begin, endelement);
    return {element, endelement};
}

auto parseNumAtoms(string::iterator begin, string::iterator end) -> pair<double, string::iterator>
{
    if(begin == end) return {1.0, begin};
    if(!isdigit(*begin)) return {1.0, begin};
    auto endnumber = std::find_if(begin, end, [](char c){return !isdigit(c);});
    double number = atof(string(begin, endnumber).c_str());
    return {number, endnumber};
}

auto findMatchedParenthesis(string::iterator begin, string::iterator end) -> string::iterator
{
    if(begin == end) return end;
    int level = 0;
    for(auto iter = begin+1; iter < end; ++iter)
    {
        level = (*iter == '(') ? level + 1 : level;
        level = (*iter == ')') ? level - 1 : level;
        if(*iter == ')' && level == -1)
            return iter;
    }
    return end;
}

auto parseFormula(string::iterator begin, string::iterator end, unordered_map<string, double>& result, double scalar) -> void
{
    if(begin == end) return;

    if(*begin != '(' && *begin != '.' && !isupper(*begin))
    {
        parseFormula(begin + 1, end, result, scalar);
    }
    else if(*begin == '(')
    {
        string::iterator begin1 = begin + 1;
        string::iterator end1   = findMatchedParenthesis(begin, end);

        auto res = parseNumAtoms(end1 + 1, end);

        const double number = res.first;

        string::iterator begin2 = res.second;
        string::iterator end2   = end;

        parseFormula(begin1, end1, result, scalar * number);
        parseFormula(begin2, end2, result, scalar);
    }
    else if(*begin == '.')
    {
        auto res = parseNumAtoms(begin + 1, end);

        const double number = res.first;

        string::iterator begin1 = res.second;
        string::iterator end1   = end;

        parseFormula(begin1, end1, result, scalar * number);
    }
    else
    {
        auto res1 = parseElementAtom(begin, end);
        auto res2 = parseNumAtoms(res1.second, end);

        string element = res1.first;
        double natoms  = res2.first;

        if(result.count(element)) result[element] += scalar * natoms;
        else result[element] = scalar * natoms;

        string::iterator begin1 = res2.second;
        string::iterator end1   = end;

        parseFormula(begin1, end1, result, scalar);
    }
}

} // namespace internal

auto elements() -> std::vector<std::string>
{
    return internal::elements;
}

auto elements(std::string formula) -> std::unordered_map<std::string, double>
{
    std::unordered_map<std::string, double> result;

    internal::parseFormula(formula.begin(), formula.end(), result, 1.0);

    return result;
}

auto atomicMass(std::string element) -> double
{
    auto entry = internal::atomicWeights.find(element);
    if(entry == internal::atomicWeights.end())
    {
        Exception exception;
        exception.error << "Cannot calculate the atomic mass of the provided chemical element.";
        exception.reason << "Element symbol `" << element << "` is unknown.";
        RaiseError(exception);
    }
    const auto gram_to_kilogram = 0.001;
    const auto atomic_mass = entry->second;
    return atomic_mass * gram_to_kilogram;
}

auto molarMass(const std::unordered_map<std::string, double>& elements) -> double
{
    double molar_mass = 0.0;
    for(const auto& pair : elements)
        molar_mass += pair.second * atomicMass(pair.first);
    return molar_mass;
}

auto molarMass(std::string formula) -> double
{
    return molarMass(elements(formula));
}

auto charge(const std::string& formula) -> double
{
    std::size_t ipos = formula.find_last_of('+');
    std::size_t ineg = formula.find_last_of('-');
    std::size_t imin = std::min(ipos, ineg);

    if(imin == std::string::npos)
        return 0.0;

    int sign = (imin == ipos) ? +1 : -1;

    if(imin + 1 == formula.size())
        return sign;

    std::string digits = formula.substr(imin + 1);

    return sign * std::stod(digits);
}

} // namespace Reaktoro
