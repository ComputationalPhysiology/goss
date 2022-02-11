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

        self._solver = solver
        self._ode = ode

    @property
    def num_nodes(self):
        return self._cpp_object.num_nodes()

    @property
    def solver(self):
        return self._solver

    @property
    def ode(self):
        return self._ode

    @property
    def num_threads(self) -> int:
        return self._cpp_object.get_num_threads()

    @num_threads.setter
    def num_threads(self, num_threads: int) -> None:
        self._cpp_object.set_num_threads(num_threads)
