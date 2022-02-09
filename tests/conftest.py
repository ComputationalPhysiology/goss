from pathlib import Path

import pytest
from goss import ODE
from gotran import load_ode

here = Path(__file__).parent.absolute()


@pytest.fixture(scope="session")
def oscilator():
    return ODE(load_ode(here.joinpath("oscilator.ode")))
