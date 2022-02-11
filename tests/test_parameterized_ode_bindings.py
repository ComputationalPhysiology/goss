from pathlib import Path

import goss
import gotran
import numpy as np
import pytest


here = Path(__file__).parent.absolute()


def test_parameterized_ode_constructors(oscilator_ode):
    ode = goss.ParameterizedODE(oscilator_ode)
    ode_copy = ode.copy()
    ode_path = goss.ParameterizedODE(here.joinpath("oscilator.ode"))
    assert ode.num_states == ode_copy.num_states == ode_path.num_states

    ode_new = goss.ParameterizedODE(num_states=3, num_parameters=1)
    assert ode_new.num_states == 3
    assert ode_new.num_parameters == 1
    assert ode.num_monitored == 0


def test_parameterized_ode_with_field_states(oscilator_ode):
    field_state_names = ["y"]
    ode = goss.ParameterizedODE(
        oscilator_ode,
        field_states=field_state_names,
    )
    assert ode.num_field_states == 1
    assert ode.num_states == 2
    assert ode.num_field_parameters == 0
    assert ode.num_parameters == 2
    assert ode.num_monitored == 0
    assert ode.state_names == ["x", "y"]
    assert ode.parameter_names == ["a", "b"]
    assert ode.field_state_names == field_state_names
    assert ode.field_state_indices == [1]
    assert ode.field_parameter_names == []
    assert ode.monitored_names == []


def test_parameterized_ode_with_invalid_field_states(oscilator_ode):
    with pytest.raises(gotran.GotranException):
        goss.ParameterizedODE(
            oscilator_ode,
            field_states=["a"],
        )


def test_parameterized_ode_with_invalid_field_parameters(oscilator_ode):
    with pytest.raises(gotran.GotranException):
        goss.ParameterizedODE(
            oscilator_ode,
            field_parameters=["x"],
        )


def test_parameterized_ode_with_field_parameters(oscilator_ode):
    field_parameter_names = ["a"]
    ode = goss.ParameterizedODE(
        oscilator_ode,
        field_parameters=field_parameter_names,
    )
    assert ode.num_field_states == 0
    assert ode.num_states == 2
    assert ode.num_field_parameters == 1
    assert ode.num_parameters == 2
    assert ode.num_monitored == 0
    assert ode.state_names == ["x", "y"]
    assert ode.parameter_names == ["a", "b"]
    assert ode.field_state_names == []
    assert ode.field_state_indices == []
    assert ode.field_parameter_names == field_parameter_names
    assert ode.monitored_names == []


def test_invalid_monitored(oscilator_ode):
    with pytest.raises(gotran.GotranException):
        goss.ParameterizedODE(
            oscilator_ode,
            monitored=["invalid_paramter"],
        )


def test_monitored(oscilator_ode):
    monitored_names = ["energy"]
    ode = goss.ParameterizedODE(
        oscilator_ode,
        monitored=monitored_names,
    )
    states = np.array([3.0, 5.0])
    monitored = ode.eval_monitored(states, 0)
    assert np.isclose(monitored, states.sum())

    assert ode.num_field_states == 0
    assert ode.num_states == 2
    assert ode.num_field_parameters == 0
    assert ode.num_parameters == 2
    assert ode.num_monitored == 1
    assert ode.state_names == ["x", "y"]
    assert ode.parameter_names == ["a", "b"]
    assert ode.field_state_names == []
    assert ode.field_state_indices == []
    assert ode.field_parameter_names == []
    assert ode.monitored_names == monitored_names
