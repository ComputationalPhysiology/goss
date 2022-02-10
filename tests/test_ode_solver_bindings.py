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


@pytest.mark.parametrize("Solver", goss.goss_non_adaptive_solvers)
def test_is_adaptive_False(Solver):
    solver = Solver()
    assert solver.is_adaptive is False


@pytest.mark.parametrize("Solver", goss.goss_adaptive_solvers)
def test_is_adaptive_True(Solver):
    solver = Solver()
    assert solver.is_adaptive is True


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_get_ode(Solver, oscilator):
    solver = Solver(oscilator)
    ode = solver.get_ode()
    assert ode.num_states == oscilator.num_states
    # FIXME: Perhaps we should assert some more
    # or implement __eq__ for ODE


@pytest.mark.parametrize("Solver", goss.goss_explicit_solvers)
def test_forward(Solver, oscilator):
    y = oscilator.get_ic()
    t = 0
    interval = 2.0
    solver = Solver(oscilator)
    y_next = y.copy()
    solver.forward(y_next, t, interval)
    assert np.allclose(y_next, [1, 2])


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_parameters(Solver):
    solver = Solver()
    parameters = solver.parameters

    # All solvers should have this paramameter
    assert np.isclose(parameters["ldt"], -1.0)

    for name, ptype in solver.parameter_names.items():

        if ptype is int:
            new_value = 42
        elif ptype is float:
            new_value = 42.0
        elif ptype is bool:
            # Use a different value
            new_value = not parameters[name]

        solver.set_parameter(name, new_value)
        new_parameters = solver.parameters
        assert np.isclose(new_parameters[name], new_value), name


@pytest.mark.parametrize("Solver", goss.goss_solvers)
def test_solve_oscilator(Solver, oscilator):
    solver = Solver(oscilator)
    u, time = solver.solve(0, 10.0, dt=0.001)

    # FIXME: Don't know why we need such a high tolerance here
    assert np.isclose(u[:, 0], np.cos(time), rtol=0.14).all()
    assert np.isclose(u[:, 1], np.sin(time), rtol=0.02).all()


@pytest.mark.parametrize("Solver", goss.goss_implicit_solvers)
def test_num_jac_comp(Solver, oscilator):
    solver = Solver(oscilator)
    assert solver.num_jac_comp() == 0


def test_AdaptiveImplicitSolver_methods(oscilator):
    solver = goss.solvers.ESDIRK23a(oscilator)

    # Just make sure notthing is failing when calling these methods
    solver.get_num_accepted()
    solver.get_num_rejected()
    solver.set_single_step_mode(True)
    solver.set_tol(atol=1e-5, rtol=1e-4)
    assert np.isclose(solver.atol, 1e-5)
    assert np.isclose(solver.rtol, 1e-4)
    solver.set_iord(42)
    assert solver.iord == 42
