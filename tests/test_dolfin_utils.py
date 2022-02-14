import goss
import pytest

try:
    import dolfin

    missing_dolfin = False
except ModuleNotFoundError:
    missing_dolfin = True

from goss.dolfinutils import DOLFINODESystemSolver


@pytest.mark.skipif(missing_dolfin, reason="Test require dolfin to be installed")
def test_DOLFINODESystemSolver_single_ODE(tentusscher_2004_fields):

    params = DOLFINODESystemSolver.default_parameters()
    mesh = dolfin.UnitCubeMesh(3, 3, 3)

    system = DOLFINODESystemSolver(mesh, odes=tentusscher_2004_fields, params=params)
    breakpoint()
