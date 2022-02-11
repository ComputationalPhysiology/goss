import goss


def test_construct_ODESysystemSolver(oscilator):
    solver = goss.solvers.ExplicitEuler()
    system_solver = goss.ODESystemSolver(1, solver, oscilator)

    ode = solver.get_ode()
    # FIXME: The ODE retured here is just an ODE and not a
    # ParameterizedODE
    assert ode.num_states == oscilator.num_states
    assert system_solver is not None
