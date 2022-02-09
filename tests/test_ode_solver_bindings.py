import goss
import numpy as np
import pytest


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_ExplicitEuler_constructor_ode(Solver, oscilator):
    solver = Solver(oscilator)
    assert solver.num_states == oscilator.num_states


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_constructor_empty(Solver):
    solver = Solver()
    assert solver.num_states == 0


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_is_adaptive(Solver):
    solver = Solver()
    assert solver.is_adaptive is False


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_get_ode(Solver, oscilator):
    solver = Solver(oscilator)
    ode = solver.get_ode()
    assert ode.num_states == oscilator.num_states
    # FIXME: Perhaps we should assert some more
    # or implement __eq__ for ODE


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_forward(Solver, oscilator):
    y = oscilator.get_ic()
    t = 0
    interval = 2.0
    solver = Solver(oscilator)
    y_next = solver.forward(y, t, interval)
    assert np.allclose(y_next, [1, 2])


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_parameters(Solver):
    solver = Solver()
    parameters = solver.parameters
    assert np.isclose(parameters["ldt"], -1.0)

    new_value = 1.0
    solver.set_parameter("ldt", new_value)
    new_parameters = solver.parameters
    assert np.isclose(new_parameters["ldt"], new_value)
