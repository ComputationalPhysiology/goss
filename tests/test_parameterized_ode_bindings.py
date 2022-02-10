from pathlib import Path

import goss
import gotran


here = Path(__file__).parent.absolute()


def test_parameterized_ode_constructors():
    ode = goss.ParameterizedODE(gotran.load_ode(here.joinpath("oscilator.ode")))
    ode_copy = ode.copy()
    ode_path = goss.ParameterizedODE(here.joinpath("oscilator.ode"))
    assert ode.num_states == ode_copy.num_states == ode_path.num_states

    ode_new = goss.ParameterizedODE(num_states=3, num_parameters=1)
    assert ode_new.num_states == 3
    assert ode_new.num_parameters == 1
