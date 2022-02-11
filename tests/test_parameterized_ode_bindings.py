from pathlib import Path

import goss
import gotran
import numpy as np
import pytest


here = Path(__file__).parent.absolute()


def test_parameterized_ode_constructors():
    ode = goss.ParameterizedODE(gotran.load_ode(here.joinpath("oscilator.ode")))
    ode_copy = ode.copy()
    ode_path = goss.ParameterizedODE(here.joinpath("oscilator.ode"))
    assert ode.num_states == ode_copy.num_states == ode_path.num_states

    ode_new = goss.ParameterizedODE(num_states=3, num_parameters=1)
    assert ode_new.num_states == 3
    assert ode_new.num_parameters == 1
    assert ode.num_monitored == 0


def test_contstruct_parameterized_ode_with_field_states():
    ode = goss.ParameterizedODE(
        gotran.load_ode(here.joinpath("oscilator.ode")),
        field_states=["x"],
    )
    assert ode.num_field_states == 1
    assert ode.num_states == 2
    assert ode.num_field_parameters == 0
    assert ode.num_parameters == 2
    assert ode.num_monitored == 0


def test_contstruct_parameterized_ode_with_invalid_field_states():
    with pytest.raises(gotran.GotranException):
        goss.ParameterizedODE(
            gotran.load_ode(here.joinpath("oscilator.ode")),
            field_states=["a"],
        )


def test_contstruct_parameterized_ode_with_invalid_field_parameters():
    with pytest.raises(gotran.GotranException):
        goss.ParameterizedODE(
            gotran.load_ode(here.joinpath("oscilator.ode")),
            field_parameters=["x"],
        )


def test_contstruct_parameterized_ode_with_field_parameters():
    ode = goss.ParameterizedODE(
        gotran.load_ode(here.joinpath("oscilator.ode")),
        field_parameters=["a"],
    )
    assert ode.num_field_states == 0
    assert ode.num_states == 2
    assert ode.num_field_parameters == 1
    assert ode.num_parameters == 2
    assert ode.num_monitored == 0


def test_invalid_monitored():
    with pytest.raises(gotran.GotranException):
        goss.ParameterizedODE(
            gotran.load_ode(here.joinpath("oscilator.ode")),
            monitored=["invalid_paramter"],
        )


def test_monitored():
    ode = goss.ParameterizedODE(
        gotran.load_ode(here.joinpath("oscilator.ode")),
        monitored=["energy"],
    )
    assert ode.num_monitored == 1
    states = np.array([3.0, 5.0])
    monitored = ode.eval_monitored(states, 0)
    assert np.isclose(monitored, states.sum())
