from pathlib import Path

import goss
import gotran
import numpy as np
import pytest


here = Path(__file__).parent.absolute()


@pytest.mark.parametrize(
    "field_states, field_state_indices",
    [([], []), (["x"], [0]), (["y"], [1]), (["x", "y"], [0, 1]), (["y", "x"], [1, 0])],
)
def test_parameterized_ode_field_states(
    field_states,
    field_state_indices,
    oscilator_ode,
):

    ode = goss.ParameterizedODE(oscilator_ode, field_states=field_states)
    assert ode.num_field_states == len(field_states)
    assert ode.field_state_names == field_states

    assert ode.num_states == 2
    assert ode.num_field_parameters == 0
    assert ode.num_parameters == 2
    assert ode.num_monitored == 0
    assert ode.state_names == ["x", "y"]
    assert ode.parameter_names == ["a", "b"]
    assert ode.field_state_indices == field_state_indices
    assert ode.field_parameter_names == []
    assert ode.monitored_names == []


@pytest.mark.parametrize("field_parameters", [[], ["a"], ["b"], ["a", "b"]])
def test_parameterized_ode_field_parameters(field_parameters, oscilator_ode):

    ode = goss.ParameterizedODE(oscilator_ode, field_parameters=field_parameters)
    assert ode.num_field_parameters == len(field_parameters)
    assert ode.field_parameter_names == field_parameters

    assert ode.num_field_states == 0
    assert ode.num_states == 2
    assert ode.num_parameters == 2
    assert ode.num_monitored == 0
    assert ode.state_names == ["x", "y"]
    assert ode.parameter_names == ["a", "b"]
    assert ode.field_state_names == []
    assert ode.field_state_indices == []
    assert ode.monitored_names == []


def test_parameterized_ode_parameters(oscilator_ode):
    ode = goss.ParameterizedODE(oscilator_ode)

    assert ode.num_parameters == oscilator_ode.num_parameters == 2

    a = ode.get_parameter("a")
    assert np.isclose(a, 1.0)
    b = ode.get_parameter("b")
    assert np.isclose(b, 1.0)

    parameters = ode.parameters
    assert np.isclose(parameters["a"], 1.0)
    assert np.isclose(parameters["b"], 1.0)

    ode.set_parameter("a", 42.0)
    ode.set_parameter("b", 13.0)

    a = ode.get_parameter("a")
    assert np.isclose(a, 42.0)
    b = ode.get_parameter("b")
    assert np.isclose(b, 13.0)

    parameters = ode.parameters
    assert np.isclose(parameters["a"], 42.0)
    assert np.isclose(parameters["b"], 13.0)


def test_parameterized_ode_constructors(oscilator_ode):
    ode = goss.ParameterizedODE(oscilator_ode)
    ode_copy = ode.copy()
    ode_path = goss.ParameterizedODE(here.joinpath("oscilator.ode"))
    assert ode.num_states == ode_copy.num_states == ode_path.num_states

    ode_new = goss.ParameterizedODE(num_states=3, num_parameters=1)
    assert ode_new.num_states == 3
    assert ode_new.num_parameters == 1
    assert ode.num_monitored == 0


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
