from pathlib import Path

import pytest
from goss import ParameterizedODE
from gotran import load_ode

here = Path(__file__).parent.absolute()


@pytest.fixture(scope="session")
def oscilator():
    return ParameterizedODE(load_ode(here.joinpath("oscilator.ode")))


@pytest.fixture(scope="session")
def tentusscher_2004():
    return ParameterizedODE(load_ode(here.joinpath("tentusscher_2004_mcell.ode ")))
