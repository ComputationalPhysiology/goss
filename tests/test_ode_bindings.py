from pathlib import Path

import goss
import gotran
import numpy as np
import pytest


here = Path(__file__).parent.absolute()


def test_ode_num_states(oscilator):
    assert oscilator.num_states == 2


def test_get_ic(oscilator):
    ic = oscilator.get_ic()
    assert np.allclose(ic, np.array([1.0, 0]))


def test_ode_eval_full(oscilator):
    values = oscilator.eval(np.array([0.1, 0.2]), 0)
    assert np.allclose(values, np.array([-0.2, 0.1]))


def test_ode_eval_index(oscilator):
    assert np.isclose(oscilator.eval_component(0, np.array([0.1, 0.2]), 0), -0.2)
    assert np.isclose(oscilator.eval_component(1, np.array([0.1, 0.2]), 0), 0.1)


def test_ode_eval_index_out_of_bounds_raises_RuntimeError(oscilator):
    with pytest.raises(RuntimeError):
        oscilator.eval_component(3, np.array([0.1, 0.2]), 0)


def test_ode_eval_negative_index_raises_TypeError(oscilator):
    with pytest.raises(TypeError):
        oscilator.eval_component(-1, np.array([0.1, 0.2]), 0)


def test_ode_compute_jacobian(oscilator):
    jac = oscilator.compute_jacobian(np.array([0.1, 0.2]), 0)
    assert np.allclose(jac, np.array([[0, -1], [1, 0]]))


def test_is_dae(oscilator):
    assert oscilator.is_dae is False


@pytest.mark.parametrize("only_linear", [True, False])
def test_linearized_rhs_only_linear_True(oscilator, only_linear):

    states = np.array([0.1, 0.2])
    rhs, linear = oscilator.linearized_eval(states, 0, only_linear)

    if only_linear:
        assert linear is None
    else:
        assert np.allclose(linear, 0.0)
    assert np.allclose(rhs, np.array([-0.2, 0.1]))


def test_ode_contructors():
    ode = goss.ODE(gotran.load_ode(here.joinpath("oscilator.ode")))
    ode_copy = ode.copy()
    ode_path = goss.ODE(here.joinpath("oscilator.ode"))
    assert ode.num_states == ode_copy.num_states == ode_path.num_states
