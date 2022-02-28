import numpy as np

from .ode import ParameterizedODE
from .solvers import ODESolver


def has_shape(y: np.ndarray, xdim: int, ydim: int) -> bool:
    if y.shape == (xdim, ydim):
        return True
    if y.size == xdim * ydim:
        return True
    return False


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
        self._tangled_storage = True

    @property
    def num_nodes(self):
        """Get the number of nodes"""
        return self._cpp_object.num_nodes()

    @property
    def solver(self) -> ODESolver:
        """Get the ODESolver"""
        return self._solver

    @property
    def ode(self) -> ParameterizedODE:
        """Get the ParameterizedODE"""
        return self._ode

    @property
    def num_threads(self) -> int:
        """Get the number of threads"""
        return self._cpp_object.get_num_threads()

    @num_threads.setter
    def num_threads(self, num_threads: int) -> None:
        """Set the number of threads"""
        self._cpp_object.set_num_threads(num_threads)

    def reset_default(self) -> None:
        """Use initial condition and initial values of
        field parameters to reset System variables
        """
        self._cpp_object.reset_default()

    @property
    def states(self) -> np.ndarray:
        """Get the whole states data array"""
        return self._cpp_object.states()

    @property
    def field_parameters(self) -> np.ndarray:
        field_parameters = np.zeros((self.num_nodes, self.ode.num_field_parameters))
        self._cpp_object.get_field_parameters(field_parameters, self._tangled_storage)
        return field_parameters

    @field_parameters.setter
    def field_parameters(self, field_parameters: np.ndarray) -> None:
        field_parameters = np.ascontiguousarray(field_parameters)
        if not has_shape(
            field_parameters,
            self.num_nodes,
            self.ode.num_field_parameters,
        ):
            raise ValueError(
                (
                    "Expected shape of field parameters to be "
                    f"{(self.num_nodes, self.ode.num_field_parameters)}, "
                    f"got {field_parameters.shape}."
                ),
            )
        self._cpp_object.set_field_parameters(field_parameters, self._tangled_storage)

    @property
    def field_states(self) -> np.ndarray:
        field_states = np.zeros((self.num_nodes, self.ode.num_field_states))
        self._cpp_object.get_field_states(field_states, self._tangled_storage)
        return field_states

    @field_states.setter
    def field_states(self, field_states: np.ndarray) -> None:
        field_states = np.ascontiguousarray(field_states)
        if not has_shape(field_states, self.num_nodes, self.ode.num_field_states):
            raise ValueError(
                (
                    "Expected shape of field states to be "
                    f"{(self.num_nodes, self.ode.num_field_states)}, "
                    f"got {field_states.shape}."
                ),
            )
        self._cpp_object.set_field_states(field_states, self._tangled_storage)

    def forward(self, t: float, interval: float) -> None:
        self._cpp_object.forward(t, float(interval))

    def solve(
        self,
        t: np.ndarray,
    ) -> np.ndarray:
        """Solve the ode for a given number of time steps

        Parameters
        ----------
        t : np.ndarray
            The time steps
        y0 : Optional[np.ndarray], optional
            Initial conditions. If not provided (default), then
            the default initial conditions will be used.

        Returns
        -------
        np.ndarray
            The states at each time point.
        """
        num_steps = t.size
        f0 = self.field_states
        field_states = np.zeros((num_steps, self.num_nodes, self.ode.num_field_states))
        field_states[0, :, :] = f0
        self._cpp_object.solve(field_states, t, num_steps, self._tangled_storage)
        return field_states
