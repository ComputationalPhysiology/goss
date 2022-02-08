import numpy as np
import pytest


def test_ode_num_states(oscilator):
    assert oscilator.num_states() == 2


def test_ode_eval_full(oscilator):

    values = np.zeros(2)
    oscilator.eval(np.array([0.1, 0.2]), 0, values)
    assert np.allclose(values, np.array([-0.2, 0.1]))


def test_ode_eval_index(oscilator):
    assert np.isclose(oscilator.eval(0, np.array([0.1, 0.2]), 0), -0.2)
    assert np.isclose(oscilator.eval(1, np.array([0.1, 0.2]), 0), 0.1)


def test_ode_eval_index_out_of_bounds_raises_RuntimeError(oscilator):
    with pytest.raises(RuntimeError):
        oscilator.eval(3, np.array([0.1, 0.2]), 0)


def test_ode_eval_negative_index_raises_TypeError(oscilator):
    with pytest.raises(TypeError):
        oscilator.eval(-1, np.array([0.1, 0.2]), 0)


def test_ode_compute_jacobian(oscilator):
    jac = np.zeros((2, 2))
    oscilator.compute_jacobian(np.array([0.1, 0.2]), 0, jac)
    assert np.allclose(jac, np.array([[0, -1], [1, 0]]))


def test_is_dae(oscilator):
    assert oscilator.is_dae() is False


def test_linearized_rhs_only_linear_True(oscilator):
    # Add some fill values to linear
    fill = 42
    linear = fill * np.ones(2)

    rhs = np.zeros(2)
    states = np.array([0.1, 0.2])
    oscilator.linearized_eval(states, 0, linear, rhs, True)

    # Linear has not changed
    assert np.allclose(linear, fill)
    assert np.allclose(rhs, np.array([-0.2, 0.1]))


def test_linearized_rhs_only_linear_False(oscilator):
    # Add some fill values to linear
    fill = 42
    linear = fill * np.ones(2)

    rhs = np.zeros(2)
    states = np.array([0.1, 0.2])
    oscilator.linearized_eval(states, 0, linear, rhs, False)

    # Now linear has changed
    assert np.allclose(linear, 0.0)
    assert np.allclose(rhs, np.array([-0.2, 0.1]))
