from core import EquilibriumSolver, EquilibriumResult, ChemicalState

_EquilibriumSolver_solve = EquilibriumSolver.solve

def _EquilibriumSolver_solve_new(solver, state, be):
    if type(state) == ChemicalState:
        return _EquilibriumSolver_solve(solver, state, be)
    else:
        result = EquilibriumResult()
        for i in xrange(len(state)):
            result += _EquilibriumSolver_solve(solver, state[i], be[i])
            if not result.succeeded:
                raise RuntimeError("Failed to calculate equilibrium state at state number %d." % i)
        return result

EquilibriumSolver.solve = _EquilibriumSolver_solve_new
