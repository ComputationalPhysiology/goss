import goss
import numpy as np
import pytest

try:
    import dolfin
    from goss import dolfinutils

    missing_dolfin = False
except ModuleNotFoundError:
    missing_dolfin = True

require_dolfin = pytest.mark.skipif(
    missing_dolfin,
    reason="Test require dolfin to be installed",
)


@pytest.fixture(scope="session")
def unitcube():
    return dolfin.UnitCubeMesh(3, 3, 3)


@pytest.fixture(scope="session")
def V(unitcube):
    return dolfin.FunctionSpace(unitcube, "Lagrange", 1)


@require_dolfin
@pytest.mark.fenics
def test_entity_to_dofs(V):
    entity_to_dof = dolfinutils.entity_to_dofs(V)
    assert len(np.unique(entity_to_dof)) == V.dim()


@require_dolfin
@pytest.mark.fenics
def test_setup_dofs_one_field(V):
    field_name = "field1"
    dofs = dolfinutils.setup_dofs(V=V, field_names=[field_name])
    assert dofs.num_dofs == {0: 64}
    assert len(dofs.dolfin_values) == 64
    assert 0 in dofs.goss_indices
    assert field_name in dofs.goss_indices[0]
    assert len(dofs.goss_indices[0][field_name]) == 64


@require_dolfin
@pytest.mark.fenics
def test_DOLFINODESystemSolver_single_ODE(tentusscher_2004_ode):

    params = dolfinutils.DOLFINODESystemSolver.default_parameters()
    mesh = dolfin.UnitCubeMesh(3, 3, 3)
    ode = goss.ParameterizedODE(
        tentusscher_2004_ode,
        field_states=["V"],
        field_parameters=["g_CaL"],
    )

    ode_solver = dolfinutils.DOLFINODESystemSolver(
        mesh,
        odes=ode,
        params=params,
    )

    V_ic = ode.get_ic()[ode.state_names.index("V")]
    solution = dolfin.Function(ode_solver.state_space)
    solution.vector()[:] = V_ic

    ode_solver.step((0, 0.1), solution)
    assert np.allclose(solution.vector().get_local(), V_ic)


@require_dolfin
@pytest.mark.fenics
def test_DOLFINODESystemSolver_muliple_ODEs(tentusscher_2004_ode, fitzhughnagumo_ode):

    params = dolfinutils.DOLFINODESystemSolver.default_parameters()
    mesh = dolfin.UnitCubeMesh(3, 3, 3)

    # Set up a left and right domain
    left = dolfin.CompiledSubDomain("x[0] <= 0.5")
    # break
    domains = dolfin.MeshFunction("size_t", mesh, 0)
    domains.set_all(0)
    left.mark(domains, 1)

    # Set up two odes
    ode0 = goss.ParameterizedODE(
        fitzhughnagumo_ode,
        field_states=["V"],
        field_parameters=["a"],
    )

    ode1 = goss.ParameterizedODE(
        tentusscher_2004_ode,
        field_states=["V"],
        field_parameters=["g_Kr", "g_Na", "g_Ks"],
    )

    odes = {0: ode0, 1: ode1}

    ode_solver = dolfinutils.DOLFINODESystemSolver(
        mesh,
        odes=odes,
        domains=domains,
        params=params,
        space="P_1",
    )

    V_ic = -85
    solution = dolfin.Function(ode_solver.state_space)
    solution.vector()[:] = V_ic  # V in mv

    ode_solver.step((0, 0.1), solution)

    # All values should be pretty close
    assert np.allclose(solution.vector().get_local(), V_ic, atol=0.1)

    # But there should be some variation
    assert np.max(np.abs(np.diff(solution.vector().get_local()))) > 1e-4
