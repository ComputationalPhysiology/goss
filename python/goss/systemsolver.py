from .ode import ParameterizedODE
from .solvers import ODESolver


class ODESystemSolver:
    def __init__(self, num_nodes: int, solver: ODESolver, ode: ParameterizedODE):

        from . import _gosscpp  # type: ignore

        self._cpp_object = _gosscpp.ODESystemSolver(
            num_nodes,
            solver._cpp_object,
            ode._cpp_object,
        )
