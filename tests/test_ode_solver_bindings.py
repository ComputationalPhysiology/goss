import goss
import pytest


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_ExplicitEuler_constructor_ode(Solver, oscilator):
    solver = Solver(oscilator)
    assert solver.num_states == oscilator.num_states


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_constructor_empty(Solver):
    solver = Solver()
    assert solver.num_states == 0
