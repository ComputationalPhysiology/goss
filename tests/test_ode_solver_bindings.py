import goss
import numpy as np
import pytest


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_solver_constructor_ode(Solver, oscilator):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls(oscilator)
    assert solver.num_states == oscilator.num_states


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_solver_constructor_ode_attach(Solver, oscilator):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    solver.attach(oscilator)
    assert solver.num_states == oscilator.num_states


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_constructor_empty(Solver):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    assert solver.num_states == 0


def test_invalid_argument_constructor():
    with pytest.raises(ValueError):
        goss.solvers.ExplicitEuler("invalid_argument")


@pytest.mark.parametrize("Solver", goss.solvers.GOSSNonAdaptiveSolvers)
def test_is_adaptive_False(Solver):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    assert solver.is_adaptive is False


@pytest.mark.parametrize("Solver", goss.solvers.GOSSAdaptiveSolvers)
def test_is_adaptive_True(Solver):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    assert solver.is_adaptive is True


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_get_ode(Solver, oscilator):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls(oscilator)
    ode = solver.get_ode()
    assert ode.num_states == oscilator.num_states
    # FIXME: Perhaps we should assert some more
    # or implement __eq__ for ODE


@pytest.mark.parametrize("Solver", goss.solvers.GOSSExplicitSolvers)
def test_forward(Solver, oscilator):
    cls = goss.solvers.solver_mapper[Solver.name]
    y = oscilator.get_ic()
    t = 0
    interval = 0.1
    solver = cls(oscilator)
    y_next = y.copy()
    solver.forward(y_next, t, interval)
    assert np.allclose(y_next, [1, 0.1], rtol=0.01)


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_set_parameter(Solver):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    parameters = solver.parameters

    for name, default_value in cls.default_parameters().items():

        if isinstance(default_value, bool):
            # Use a different value
            new_value = not parameters[name]
        elif isinstance(default_value, int):
            new_value = 42
        elif isinstance(default_value, float):
            new_value = 42.0

        solver.set_parameter(name, new_value)
        new_parameters = solver.parameters
        assert np.isclose(new_parameters[name], new_value), name


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_update_parameters(Solver):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    parameters = solver.parameters

    new_parameters = {}
    for name, default_value in cls.default_parameters().items():

        if isinstance(default_value, bool):
            # Use a different value
            new_value = not parameters[name]
        elif isinstance(default_value, int):
            new_value = 42
        elif isinstance(default_value, float):
            new_value = 42.0

        new_parameters[name] = new_value

    solver.update_parameters(new_parameters)
    for name, new_value in new_parameters.items():
        value = solver.get_parameter(name)
        assert np.isclose(value, new_value), name


def test_set_invalid_parameter_type():
    solver = goss.solvers.ThetaSolver()
    with pytest.raises(TypeError):
        solver.set_parameter("theta", "hello")


def test_set_invalid_parameter_name():
    solver = goss.solvers.ThetaSolver()
    with pytest.raises(KeyError):
        solver.set_parameter("invalid_parmeter", 42.0)


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_solve_oscilator(Solver, oscilator):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls(oscilator)
    solver.internal_time_step = 0.0001
    dt = 0.1
    t = np.arange(0, 10.0 + dt, dt)
    u = solver.solve(t)

    # FIXME: Don't know why we need such a high tolerance here
    assert np.isclose(u[:, 0], np.cos(t), rtol=0.14).all()
    assert np.isclose(u[:, 1], np.sin(t), rtol=0.02).all()


@pytest.mark.parametrize("Solver", goss.solvers.GOSSImplicitSolvers)
def test_num_jac_comp(Solver, oscilator):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls(oscilator)
    assert solver.num_jac_comp() == 0


def test_AdaptiveImplicitSolver_methods(oscilator):
    solver = goss.solvers.ESDIRK23a(oscilator)

    # Just make sure notthing is failing when calling these methods
    solver.get_num_accepted()
    solver.get_num_rejected()
    assert np.isclose(solver.get_current_time(), 0.0)
    assert np.isclose(solver.get_current_time_step(), 0.0)
    solver.set_single_step_mode(True)
    solver.set_tol(atol=1e-5, rtol=1e-4)
    assert np.isclose(solver.atol, 1e-5)
    assert np.isclose(solver.rtol, 1e-4)
    solver.set_iord(42)
    assert solver.iord == 42


@pytest.mark.parametrize("Solver", goss.solvers.GOSSAdaptiveSolvers)
def test_Adaptive_methods(Solver):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    solver.reset()
    assert solver.nfevals == 0
    assert solver.ndtsa == 0
    assert solver.ndtsr == 0


def test_repr_str():
    solver = goss.solvers.ESDIRK23a()
    assert str(solver) == "ESDIRK23a"
    assert repr(solver) == "ESDIRK23a()"


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_internal_time_step(Solver):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()
    value = 0.42
    solver.internal_time_step = value
    assert np.isclose(solver.internal_time_step, value)


@pytest.mark.parametrize("Solver", goss.solvers.GOSSImplicitSolvers)
def test_ImplicitODESolver_compute_factorized_jacobian(Solver, fitzhughnagumo_ode):
    ode = goss.ParameterizedODE(fitzhughnagumo_ode)
    cls = goss.solvers.solver_mapper[Solver.name]

    solver = cls(ode)
    assert solver.num_jac_comp() == 0
    solver.compute_factorized_jacobian(y=ode.get_ic(), t=0.0, dt=0.1, alpha=0.5)
    assert solver.num_jac_comp() == 1
