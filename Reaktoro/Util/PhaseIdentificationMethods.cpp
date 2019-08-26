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

#include "PhaseIdentificationMethods.hpp"

#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Math/Roots.hpp>

#include <Reaktoro/deps/eigen3/Eigen/Dense>
#include <Reaktoro/deps/eigen3/unsupported/Eigen/Polynomials>


namespace Reaktoro {
namespace PhaseID {

auto volumeMethod(const ThermoScalar& Temperature, const ThermoScalar& Pressure, const ChemicalScalar& Z, ChemicalScalar& b) -> PhaseType
{
    auto Volume = Z * universalGasConstant * Temperature / Pressure;
    if ((Volume.val / b.val) > 1.75)
        return PhaseType::Gas;

    return PhaseType::Liquid;
}

auto isothermalCompressibilityMethod(const ThermoScalar& Temperature, const ThermoScalar& Pressure, const ChemicalScalar& Z) -> PhaseType
{
    auto Volume = Z * universalGasConstant*Temperature / Pressure;
    auto dkdt = (1.0 / (Volume.val*Volume.val))*Volume.ddP*Volume.ddT;
    
    if (dkdt <= 0.0)
        return PhaseType::Gas;

    return PhaseType::Liquid;
}

auto pressureComparison(const ThermoScalar& Pressure, const ThermoScalar& Temperature, const ChemicalScalar& amix,
                        const ChemicalScalar& bmix, const ChemicalScalar& A, const ChemicalScalar& B, const ChemicalScalar& C,
                        const double epsilon, const double sigma) -> PhaseType
{
    
    auto p = [&](double V) -> double
    {
        return ((universalGasConstant*Temperature.val) / (V - bmix.val)) - (amix.val / ((V + epsilon * bmix.val) * (V + sigma * bmix.val)));
    };
    
    auto k1 = epsilon * bmix.val;
    auto k2 = sigma * bmix.val;

    // Computing parameters AP, BP, CP, DP and EP of equation AP*P^4 + BP*P^3 + CP*P^2 + DP*P + EP = 0,
    // which gives the values of P where the EoS changes slope (Local max and min)
    const auto R = universalGasConstant;
    const auto T = Temperature.val;
    const auto AP = R * T;
    const auto BP = 2 * R * T * (k2 + k1) - 2 * amix.val;
    const auto CP = R * T * (k2 * k2 + 4.0 * k1 * k2 + k1 * k1) - amix.val * (k1 + k2 - 4 * bmix.val);
    const auto DP = 2 * R * T * (k1 * k2 * k2 + k1 * k1 * k2) - 2 * amix.val * (bmix.val * bmix.val - k2 * bmix.val - k1 * bmix.val);
    const auto EP = R * T * k1 * k1 * k2 * k2 - amix.val * (k1 + k2) * bmix.val * bmix.val;

    auto polynomial_solver = Eigen::PolynomialSolver<double, 4>(
        Eigen::Matrix<double, 5, 1>{EP, DP, CP, BP, AP});

    constexpr auto abs_imaginary_threshold = 1e-15;
    auto real_roots = std::vector<double>();
    real_roots.reserve(5);
    polynomial_solver.realRoots(real_roots, abs_imaginary_threshold);

    // removing roots lower than bmix, no physical meaning
    auto new_end = std::remove_if(real_roots.begin(), real_roots.end(), [&](double r){
        return r < bmix.val;
    });
    real_roots.resize(new_end - real_roots.begin());
    
    if (real_roots.size() == 0) 
    {
        return PhaseType::Gas;
    }

    std::vector<double> pressures;
    for (const auto& volume : real_roots)
        pressures.push_back(p(volume));

    auto Pmin = std::min_element(pressures.begin(), pressures.end());
    auto Pmax = std::max_element(pressures.begin(), pressures.end());
    
    if (Pressure.val < *Pmin)
        return PhaseType::Gas;
    if (Pressure.val > *Pmax)
        return PhaseType::Liquid;

    Exception exception;
    exception.error << "Could not define phase type.";
    exception.reason << "gibbsEnergyAndEquationOfStateMethod has received one Z but the pressure is between Pmin and Pmax.";
    RaiseError(exception);
}


auto gibbsResidualEnergyComparison(const ThermoScalar& Pressure, const ThermoScalar& Temperature, const ChemicalScalar& amix,
                                   const ChemicalScalar& bmix, const ChemicalScalar& A, const ChemicalScalar& B,
                                   std::vector<ChemicalScalar> Zs, const double epsilon, const double sigma) -> PhaseType
{
    // Computing the values of residual Gibbs energy for all Zs
    std::vector<ChemicalScalar> Gs;
    for (const auto Z : Zs)
    {
        const double factor = -1.0 / (3 * Z.val*Z.val + 2 * A.val*Z.val + B.val);
        const ChemicalScalar beta = Pressure * bmix / (universalGasConstant * Temperature);
        const ChemicalScalar q = amix / (bmix * universalGasConstant * Temperature);
        
        // Calculate the integration factor I and its temperature derivative IT
        ChemicalScalar I;
        if (epsilon != sigma) 
            I = log((Z + sigma * beta) / (Z + epsilon * beta)) / (sigma - epsilon);
        else 
            I = beta / (Z + epsilon * beta);        

        Gs.push_back(universalGasConstant*Temperature*(Z - 1 - log(Z - beta) - q * I));
    }
    
    if (Gs[0].val < Gs[1].val)
        return PhaseType::Gas;
    else
       return PhaseType::Liquid;
    
}


auto gibbsEnergyAndEquationOfStateMethod(const ThermoScalar& Pressure, const ThermoScalar& Temperature, const ChemicalScalar& amix,
                                         const ChemicalScalar& bmix, const ChemicalScalar& A, const ChemicalScalar& B, const ChemicalScalar& C,
                                         std::vector<ChemicalScalar> Zs, const double epsilon, const double sigma) -> PhaseType
{
    if (Zs.size() == 1)
    {
        return pressureComparison(Pressure, Temperature, amix, bmix, A, B, C, epsilon, sigma);
    }
    else
    {
        return gibbsResidualEnergyComparison(Pressure, Temperature, amix, bmix, A, B, Zs, epsilon, sigma);
    }
}

}

}