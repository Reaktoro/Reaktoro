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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    ChemicalEditor editor;
    editor.addAqueousPhase("CO2(aq) CO3-- Ca(HCO3)+ Ca++ CaCO3(aq) CaCl+ CaCl2(aq) CaOH+ Cl- H+ H2O(l) HCO3- HCl(aq) Mg(HCO3)+ Mg++ MgCO3(aq) MgCl+ MgOH+ Na+ NaCl(aq) NaOH(aq) OH-");
    editor.addGaseousPhase("H2O(g) CO2(g)");

    ChemicalSystem system(editor);

    EquilibriumProblem problem0(system);
    problem0.setTemperature(200, "celsius");
    problem0.setPressure(300, "bar");
    problem0.add("H2O", 1, "kg");
    problem0.add("CO2", 1, "mol");
    problem0.add("NaCl", 1, "mol");
    problem0.add("CaCl2", 0.1, "mol");
    problem0.add("MgCl2", 0.05, "mol");

    const auto state0 = equilibrate(problem0);
    const auto props0 = state0.properties();

    const auto U = props0.internalEnergy().val;
    const auto V = props0.volume().val;

    EquilibriumInverseProblem problem(system);
    problem.setTemperature(180.0, "celsius");
    problem.setPressure(270.0, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 1, "mol");
    problem.add("NaCl", 1, "mol");
    problem.add("CaCl2", 0.1, "mol");
    problem.add("MgCl2", 0.05, "mol");

    problem.unknownTemperature();
    problem.unknownPressure();

    problem.fixVolume(V, "m3");
    problem.fixInternalEnergy(U, "J");

    EquilibriumOptions options;
    // options.optimum.output.active = false
    options.optimum.tolerance = 1000;
    options.nonlinear.linesearch = true;
    options.hessian = GibbsHessian::Exact;
    options.epsilon = 1e-14;
    options.nonlinear.output.active = true;

    ChemicalState state(system);

    EquilibriumInverseSolver solver(system);
    solver.setOptions(options);
    solver.solve(state, problem);

    state.output("demo-equilibrium-fixed-UV.txt");
}
