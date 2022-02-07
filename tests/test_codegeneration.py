from pathlib import Path
import pytest
from gotran import load_ode
from goss import jit

_here = Path(__file__).absolute().parent


@pytest.fixture
def ode():
    _ode = load_ode(_here.joinpath("oscilator.ode").as_posix())
    return _ode


def test_codegeneration(ode):
    module = jit(ode)