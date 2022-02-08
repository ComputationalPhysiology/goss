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
