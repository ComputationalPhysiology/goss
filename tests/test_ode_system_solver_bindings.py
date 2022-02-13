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


def test_get_field_state_components_raise_RuntimeError_when_no_fields_states(oscilator):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 3
    system_solver = goss.ODESystemSolver(num_nodes, solver, oscilator)
    with pytest.raises(RuntimeError):
        system_solver.get_field_state_components([0])


def test_get_field_state_components_multiple(tentusscher_2004_fields):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 3
    system_solver = goss.ODESystemSolver(num_nodes, solver, tentusscher_2004_fields)

    ic = tentusscher_2004_fields.get_ic()
    field_state_indices = tentusscher_2004_fields.field_state_indices
    assert len(field_state_indices) == 2
    field_state_values = ic[field_state_indices]

    states = system_solver.states
    assert states.shape == (num_nodes, tentusscher_2004_fields.num_states)

    field_states_comp = system_solver.get_field_state_components([0, 1])
    assert (field_states_comp[:, 0] == field_state_values[0]).all()
    assert (field_states_comp[:, 1] == field_state_values[1]).all()

    field_states_comp = system_solver.get_field_state_components([1, 0])
    assert (field_states_comp[:, 0] == field_state_values[1]).all()
    assert (field_states_comp[:, 1] == field_state_values[0]).all()


@pytest.mark.parametrize(
    "component",
    [
        0,
        pytest.param(1, marks=pytest.mark.xfail(reason="some bug")),
    ],
)
def test_get_field_state_components_single(component, tentusscher_2004_fields):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 3
    system_solver = goss.ODESystemSolver(num_nodes, solver, tentusscher_2004_fields)

    ic = tentusscher_2004_fields.get_ic()
    field_state_indices = tentusscher_2004_fields.field_state_indices
    assert len(field_state_indices) == 2
    field_state_values = ic[field_state_indices]

    states = system_solver.states
    assert states.shape == (num_nodes, tentusscher_2004_fields.num_states)

    field_states_comp = system_solver.get_field_state_components([component])
    assert (field_states_comp == field_state_values[component]).all()


def test_field_states(tentusscher_2004_fields):
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


def test_field_parameters(tentusscher_2004_fields):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 3
    system_solver = goss.ODESystemSolver(num_nodes, solver, tentusscher_2004_fields)

    field_parameter_names = tentusscher_2004_fields.field_parameter_names
    field_parameters = system_solver.field_parameters
    parameters = tentusscher_2004_fields.parameters

    # Check default values
    for i, name in enumerate(field_parameter_names):
        assert np.allclose(field_parameters[:, i], parameters[name])

    # Test setter
    field_parameters[:] = 0
    system_solver.field_parameters = field_parameters
    zero_parameters = system_solver.field_parameters
    assert (zero_parameters == 0).all()
