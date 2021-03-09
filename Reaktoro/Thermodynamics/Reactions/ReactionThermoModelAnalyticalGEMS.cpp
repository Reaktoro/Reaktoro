// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "ReactionThermoModelAnalyticalGEMS.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

auto ReactionThermoModelAnalyticalGEMS(Param A0, Param A1, Param A2, Param A3, Param A4, Param A5, Param A6) -> ReactionThermoModel
{
    auto evalfn = [=](ReactionThermoProps& props, ReactionThermoArgs args)
    {
        ReactionThermoArgsDecl(args);
        const auto R = universalGasConstant;
        const auto T2 = T*T;
        const auto T3 = T*T2;
        const auto T05 = sqrt(T);

        props.dG0 = -R*T * (A0 + A1*T + A2/T + A3*log(T) + A4/T2 + A5*T2 + A6/T05) * ln10;
        props.dH0 = R * (A1*T2 - A2 + A3*T - 2*A4/T + 2*A5*T3 - 0.5*A6*T05) * ln10;
    };

    Params params = { A0, A1, A2, A3, A4, A5, A6 };

    return ReactionThermoModel(evalfn, params);
}

} // namespace Reaktoro