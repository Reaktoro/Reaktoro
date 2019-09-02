#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database database("supcrt98.xml");

    ChemicalEditor editor(database);
    editor.addAqueousPhaseWithElements("H O C Na Cl");
    editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
    editor.addMineralPhase("Halite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.setTemperature(60, "celsius");
    problem.setPressure(300, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 100, "g");
    problem.add("NaCl", 0.1, "mol");

    ChemicalState state = equilibrate(problem);

    state.output("state-cmake-example.txt");
}
