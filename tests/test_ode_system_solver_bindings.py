import goss
import numpy as np
import pytest


@pytest.fixture
def odesystemsolver(oscilator):
    solver = goss.solvers.ExplicitEuler()
    return goss.ODESystemSolver(1, solver, oscilator)


def test_construct_ODESysystemSolver(oscilator):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 1
    system_solver = goss.ODESystemSolver(num_nodes, solver, oscilator)

    assert system_solver.num_nodes == num_nodes

    ode = solver.get_ode()
    # FIXME: The ODE retured here is just an ODE and not a
    # ParameterizedODE
    assert ode.num_states == oscilator.num_states

    # Make sure we have attached the correct solver
    assert id(solver._cpp_object) == id(system_solver._cpp_object.solver())
    assert id(oscilator._cpp_object) == id(system_solver._cpp_object.ode())


def test_has_openmp():
    assert goss.has_openmp() is True


def test_num_threads(odesystemsolver):
    assert odesystemsolver.num_threads == 0

    num_threads = 1
    odesystemsolver.num_threads = num_threads
    assert odesystemsolver.num_threads == num_threads


def test_reset_default(odesystemsolver):
    # TODO: Do something
    odesystemsolver.reset_default()
    # TODO: Make sure the settings are restored


def test_states(oscilator):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 3
    system_solver = goss.ODESystemSolver(num_nodes, solver, oscilator)
    states = system_solver.states
    assert states.shape == (num_nodes, oscilator.num_states)
    ic = oscilator.get_ic()
    for node in range(num_nodes):
        assert np.allclose(states[node, :], ic)


def test_field_states(tentusscher_2004_fields: goss.ParameterizedODE):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 3
    system_solver = goss.ODESystemSolver(num_nodes, solver, tentusscher_2004_fields)

    ic = tentusscher_2004_fields.get_ic()
    field_state_indices = tentusscher_2004_fields.field_state_indices
    field_state_values = ic[field_state_indices]

    # Test getter
    field_states = system_solver.field_states
    assert field_states.shape == (num_nodes, tentusscher_2004_fields.num_field_states)
    assert (field_states[:, 0] == field_state_values[0]).all()
    assert (field_states[:, 1] == field_state_values[1]).all()

    # Test setter
    field_states[:] = 0
    system_solver.field_states = field_states
    zero_field_states = system_solver.field_states
    assert (zero_field_states == 0).all()


def test_field_parameters(tentusscher_2004_fields: goss.ParameterizedODE):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 4
    system_solver = goss.ODESystemSolver(num_nodes, solver, tentusscher_2004_fields)

    field_parameter_names = tentusscher_2004_fields.field_parameter_names
    field_parameters = system_solver.field_parameters
    parameters = tentusscher_2004_fields.parameters

    assert field_parameters.shape == (
        num_nodes,
        tentusscher_2004_fields.num_field_parameters,
    )

    # Check default values
    for i, name in enumerate(field_parameter_names):
        assert np.allclose(field_parameters[:, i], parameters[name])

    # Test setter
    field_parameters[:] = 0
    system_solver.field_parameters = field_parameters
    zero_parameters = system_solver.field_parameters
    assert (zero_parameters == 0).all()


@pytest.mark.parametrize("Solver", goss.solvers.GOSSSolvers)
def test_solve_system_oscilator(Solver, oscilator_ode):
    cls = goss.solvers.solver_mapper[Solver.name]
    solver = cls()

    ode = goss.ParameterizedODE(
        oscilator_ode,
        field_states=["x", "y"],
        field_parameters=["a", "b"],
    )

    num_nodes = 4
    system_solver = goss.ODESystemSolver(num_nodes, solver, ode)

    solver.internal_time_step = 0.0001
    dt = 0.1
    t = np.arange(0, 10.0 + dt, dt)
    field_states = system_solver.solve(t)

    # FIXME: Don't know why we need such a high tolerance here
    for i in range(num_nodes):
        assert np.isclose(field_states[:, i, 0], np.cos(t), rtol=0.14).all()
        assert np.isclose(field_states[:, i, 1], np.sin(t), rtol=0.02).all()
