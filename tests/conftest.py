from pathlib import Path
import pytest
from goss import jit
from gotran import load_ode

here = Path(__file__).parent.absolute()


@pytest.fixture(scope="session")
def oscilator():
    return jit(load_ode(here.joinpath("oscilator.ode")))
