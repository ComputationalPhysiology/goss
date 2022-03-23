import itertools as it

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
@pytest.mark.parametrize("family", ["CG", "DG", "DG_2"])
def test_entity_to_dofs(family, unitcube):
    if family == "DG_2":
        W = dolfin.VectorFunctionSpace(unitcube, family.split("_")[0], 1, dim=2)
    else:
        W = dolfin.FunctionSpace(unitcube, family, 1)
    entity_to_dof = dolfinutils.entity_to_dofs(W)
    assert len(np.unique(entity_to_dof)) == W.dim()


@require_dolfin
@pytest.mark.fenics
@pytest.mark.parametrize(
    "field_names, family",
    it.product([["field1"], ["field1", "field2"]], ["CG", "DG"]),
)
def test_setup_dofs_one_field(field_names, family, unitcube):
    degree = 1 if family == "CG" else 0
    if len(field_names) == 1:
        W = dolfin.FunctionSpace(unitcube, family, degree)
    else:
        W = dolfin.VectorFunctionSpace(unitcube, family, degree, dim=len(field_names))
    value = 42
    dim = 0 if family == "CG" else 3
    domains = dolfin.MeshFunction("size_t", unitcube, dim, value=value)
    dofs = dolfinutils.setup_dofs(V=W, domains=domains, field_names=field_names)
    num_dofs = unitcube.num_vertices() if family == "CG" else unitcube.num_cells()

    assert dofs.num_dofs == {value: num_dofs}
    assert len(dofs.dolfin_values) == W.dim()
    assert value in dofs.goss_indices
    assert field_names[0] in dofs.goss_indices[value]
    assert len(dofs.goss_indices[value][field_names[0]]) == num_dofs


@require_dolfin
@pytest.mark.fenics
def test_DOLFINParameterizedODE(tentusscher_2004_ode, V):
    fp_name = "g_Kr"
    goss_ode = goss.ParameterizedODE(
        tentusscher_2004_ode,
        field_states=["V"],
        field_parameters=[fp_name],
    )

    ode = dolfinutils.DOLFINParameterizedODE.from_ode(goss_ode)
    assert ode.num_states() == goss_ode.num_states
    # Try to set a field parameter
    assert len(ode.field_params) == 0
    assert len(ode.changed_field_parameters) == 0

    f = dolfin.Function(V)
    ode.set_parameter(fp_name, f)
    assert len(ode.field_params) == 1
    assert fp_name in ode.field_params
    assert fp_name in ode.changed_field_parameters

    ode.set_parameter("g_Ks", 13)
    assert np.isclose(ode.get_parameter("g_Ks"), 13)
    assert len(ode.field_params) == 1

    # Try to intialize the solver
    dolfinutils.DOLFINODESystemSolver(
        V.mesh(),
        odes=ode,
    )


@require_dolfin
@pytest.mark.fenics
def test_DOLFINODESystemSolver_single_ODE(tentusscher_2004_ode):

    params = dolfinutils.DOLFINODESystemSolver.default_parameters()
    mesh = dolfin.UnitCubeMesh(3, 3, 3)
    ode = dolfinutils.DOLFINParameterizedODE(
        tentusscher_2004_ode,
        field_states=["V", "Ca_i"],
        field_parameters=["g_CaL"],
    )

    V_g_CaL = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    g_CaL = dolfin.Function(V_g_CaL)
    g_CaL.assign(dolfin.Constant(ode.get_parameter("g_CaL")))
    ode.set_parameter("g_CaL", g_CaL)

    ode_solver = dolfinutils.DOLFINODESystemSolver(
        mesh,
        odes=ode,
        params=params,
    )
    ode.set_parameter("g_CaL", g_CaL)

    V_ic = ode.get_ic()[ode.state_names.index("V")]
    Cai_ic = ode.get_ic()[ode.state_names.index("Ca_i")]
    solution = ode_solver.vs

    VS0 = ode_solver.state_space.sub(0)
    VS1 = ode_solver.state_space.sub(1)
    W = VS0.collapse()
    v0_assigner = dolfin.FunctionAssigner(VS0, W)
    v1_assigner = dolfin.FunctionAssigner(VS1, W)

    V = dolfin.Function(W)
    V.assign(dolfin.Constant(V_ic))

    Cai = dolfin.Function(W)
    Cai.assign(dolfin.Constant(Cai_ic))

    v0_assigner.assign(solution.sub(0), V)
    v1_assigner.assign(solution.sub(1), Cai)

    v0_assigner_rev = dolfin.FunctionAssigner(W, VS0)
    v1_assigner_rev = dolfin.FunctionAssigner(W, VS1)

    ode_solver.step((0, 0.1))

    v0_assigner_rev.assign(V, solution.sub(0))
    v1_assigner_rev.assign(Cai, solution.sub(1))

    # Solution should be fairly close but not identical
    assert np.linalg.norm(V.vector().get_local() - V_ic) > 0
    assert np.allclose(V.vector().get_local(), V_ic, rtol=0.001)
    assert np.linalg.norm(Cai.vector().get_local() - Cai_ic) > 0
    assert np.allclose(Cai.vector().get_local(), Cai_ic, rtol=0.001)


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
    ode0 = goss.dolfinutils.DOLFINParameterizedODE(
        fitzhughnagumo_ode,
        field_states=["V"],
        field_parameters=["a"],
    )

    ode1 = goss.dolfinutils.DOLFINParameterizedODE(
        tentusscher_2004_ode,
        field_states=["V"],
        field_parameters=["g_Kr", "g_Na", "g_Ks"],
    )

    odes = {0: ode0, 1: ode1}

    V_ic = -86.2

    for ode in odes.values():
        ode.set_initial_conditions(V=V_ic)

    ode_solver = dolfinutils.DOLFINODESystemSolver(
        mesh,
        odes=odes,
        domains=domains,
        params=params,
        space="P_1",
    )

    ode_solver.step((0, 0.1))

    solution = ode_solver.vs

    # All values should be pretty close
    assert np.allclose(solution.vector().get_local(), V_ic, atol=2.0)

    # But there should be some variation
    assert np.max(np.abs(np.diff(solution.vector().get_local()))) > 1e-4
