from pathlib import Path

import goss
import numpy as np
import pytest

# import scipy.linalg


here = Path(__file__).parent.absolute()


def test_ode_num_states(oscilator):
    assert oscilator.num_states == 2


def test_get_ic(oscilator):
    ic = oscilator.get_ic()
    assert np.allclose(ic, np.array([1.0, 0]))


def test_ode_eval_full(oscilator):

    values = oscilator.eval(np.array([0.1, 0.2]), 0)
    assert np.allclose(values, np.array([-0.2, 0.1]))


def test_forward_backward_subst(tentusscher_2004):
    ic = tentusscher_2004.get_ic()
    jac = tentusscher_2004.compute_jacobian(ic, 0)
    x = np.zeros_like(ic)
    tentusscher_2004.forward_backward_subst(jac, ic, x)
    # TODO: Add some asserts here


def test_ode_eval_full_list(oscilator):
    values = oscilator.eval([0.1, 0.2], 0)
    assert np.allclose(values, np.array([-0.2, 0.1]))


def test_ode_eval_full_invalid_states_raises_ValueError(oscilator):
    with pytest.raises(ValueError):
        oscilator.eval(np.array([0.1, 0.2, 0.3]), 0)


def test_ode_eval_index(oscilator):
    assert np.isclose(oscilator.eval_component(0, np.array([0.1, 0.2]), 0), -0.2)
    assert np.isclose(oscilator.eval_component(1, np.array([0.1, 0.2]), 0), 0.1)


def test_ode_eval_index_out_of_bounds_raises_RuntimeError(oscilator):
    with pytest.raises(RuntimeError):
        oscilator.eval_component(3, np.array([0.1, 0.2]), 0)


def test_ode_eval_negative_index_raises_TypeError(oscilator):
    with pytest.raises(TypeError):
        oscilator.eval_component(-1, np.array([0.1, 0.2]), 0)


def test_ode_compute_jacobian_oscilator(oscilator):
    jac = oscilator.compute_jacobian(np.array([0.1, 0.2]), 0)
    assert np.allclose(jac, np.array([[0, -1], [1, 0]]))


def test_ode_jacobian_tentusscher_2004(tentusscher_2004):

    ic = tentusscher_2004.get_ic()
    jac = tentusscher_2004.compute_jacobian(ic, 0)
    # jac_orig = jac.copy()
    tentusscher_2004.lu_factorize(jac)
    # TODO: Any way we can verify the implementation with this?
    # p, l, u = scipy.linalg.lu(jac_orig)
    # pl, u = scipy.linalg.lu(jac_orig, permute_l=True)

    # assert np.allclose(jac, np.array([[0, -1], [1, 0]]))


def test_is_dae(oscilator):
    assert oscilator.is_dae is False


def test_monitored_values(fitzhughnagumo_ode):
    ode = goss.ParameterizedODE(fitzhughnagumo_ode, monitored=["p_42", "p_V"])
    ic = ode.get_ic()
    num_time_steps = 3
    states = np.zeros((num_time_steps, ode.num_states))
    for i in range(num_time_steps):
        states[i, :] = ic
    time = np.arange(num_time_steps)
    m = ode.monitored_values(states, time)

    assert m.shape == (num_time_steps, 2)
    assert np.allclose(m[:, 0], 42)
    assert np.allclose(m[:, 1], ic[0])


@pytest.mark.parametrize("only_linear", [True, False])
def test_linearized_rhs_only_linear_True(oscilator, only_linear):

    states = np.array([0.1, 0.2])
    rhs, linear = oscilator.linearized_eval(states, 0, only_linear)

    if only_linear:
        assert linear is None
    else:
        assert np.allclose(linear, 0.0)
    assert np.allclose(rhs, np.array([-0.2, 0.1]))


def test_ode_contructors(oscilator_ode):
    ode = goss.ODE(oscilator_ode)
    ode_copy = ode.copy()
    ode_path = goss.ODE(here.joinpath("oscilator.ode"))
    assert ode.num_states == ode_copy.num_states == ode_path.num_states


def test_invalid_ode_construtor():
    with pytest.raises(ValueError):
        goss.ODE(42)
