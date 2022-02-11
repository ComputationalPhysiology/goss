import goss


def test_construct_ODESysystemSolver(oscilator):
    solver = goss.solvers.ExplicitEuler()
    num_nodes = 1
    system_solver = goss.ODESystemSolver(num_nodes, solver, oscilator)

    assert system_solver.num_nodes == num_nodes

    ode = solver.get_ode()
    # FIXME: The ODE retured here is just an ODE and not a
    # ParameterizedODE
    assert ode.num_states == oscilator.num_states

    # Make sure we have attached the correct solver
    assert id(solver._cpp_object) == id(system_solver._cpp_object.solver())
    assert id(oscilator._cpp_object) == id(system_solver._cpp_object.ode())


def test_has_openmp():
    assert goss.has_openmp() is True


def test_num_threads(oscilator):
    solver = goss.solvers.ExplicitEuler()
    system_solver = goss.ODESystemSolver(1, solver, oscilator)
    assert system_solver.num_threads == 0

    num_threads = 1
    system_solver.num_threads = num_threads
    assert system_solver.num_threads == num_threads
