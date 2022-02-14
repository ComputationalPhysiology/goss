import pytest

try:
    import dolfin
    from goss import dolfinutils

    missing_dolfin = False
except ModuleNotFoundError:
    missing_dolfin = True


@pytest.fixture(scope="session")
def unitcube():
    return dolfin.UnitCubeMesh(3, 3, 3)


@pytest.fixture(scope="session")
def V(unitcube):
    return dolfin.FunctionSpace(unitcube, "Lagrange", 1)


@pytest.mark.skipif(missing_dolfin, reason="Test require dolfin to be installed")
def test_setup_dofs(V):
    dofs = dolfinutils.setup_dofs(V=V, field_names=["field1"])
    # assert something


@pytest.mark.skipif(missing_dolfin, reason="Test require dolfin to be installed")
def test_DOLFINODESystemSolver_single_ODE(tentusscher_2004_fields):

    params = dolfinutils.DOLFINODESystemSolver.default_parameters()
    mesh = dolfin.UnitCubeMesh(3, 3, 3)

    system = dolfinutils.DOLFINODESystemSolver(
        mesh,
        odes=tentusscher_2004_fields,
        params=params,
    )
