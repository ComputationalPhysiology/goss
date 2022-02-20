from pathlib import Path

import goss
import pytest
from gotran import load_ode

here = Path(__file__).parent.absolute()


@pytest.fixture(scope="session")
def oscilator_ode():
    return load_ode(here.joinpath("oscilator.ode"))


@pytest.fixture(scope="session")
def fitzhughnagumo_ode():
    return load_ode(here.joinpath("fitzhughnagumo.ode"))


@pytest.fixture(scope="session")
def oscilator(oscilator_ode):
    return goss.ParameterizedODE(oscilator_ode)


@pytest.fixture(scope="session")
def tentusscher_2004_ode():
    return load_ode(here.joinpath("tentusscher_2004_mcell.ode "))


@pytest.fixture(scope="session")
def tentusscher_2004(tentusscher_2004_ode):
    return goss.ParameterizedODE(tentusscher_2004_ode)


@pytest.fixture(scope="session")
def tentusscher_2004_fields(tentusscher_2004_ode):
    return goss.ParameterizedODE(
        tentusscher_2004_ode,
        field_states=["Ca_i", "V"],
        field_parameters=["g_Kr", "g_Na", "g_Ks"],
    )


@pytest.fixture(scope="session")
def tentusscher_panfilov_2006_M_cell_ode():
    return load_ode(here.joinpath("tentusscher_panfilov_2006_M_cell.ode "))
