import numpy as np

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

    def _check_field_components(self, components: np.ndarray):
        if components.max() >= self.ode.num_field_states:
            raise RuntimeError(
                (
                    f"Number of field states are {self.ode.num_field_states}. "
                    f"Cannot evalute field component number {components.max()}"
                ),
            )

    @property
    def field_parameters(self) -> np.ndarray:
        field_parameters = np.zeros((self.num_nodes, self.ode.num_field_parameters))
        self._cpp_object.get_field_parameters(field_parameters, self._tangled_storage)
        return field_parameters

    @field_parameters.setter
    def field_parameters(self, field_parameters: np.ndarray) -> None:
        field_parameters = np.ascontiguousarray(field_parameters)
        if field_parameters.shape != (self.num_nodes, self.ode.num_field_parameters):
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
        if field_states.shape != (self.num_nodes, self.ode.num_field_states):
            raise ValueError(
                (
                    "Expected shape of field states to be "
                    f"{(self.num_nodes, self.ode.num_field_states)}, "
                    f"got {field_states.shape}."
                ),
            )
        self._cpp_object.set_field_states(field_states, self._tangled_storage)

    def forward(self, t: float, interval: float) -> None:
        self._cpp_object.forward(t, interval)

    # def set_field_state_components(
    #     self,
    #     component_field_states: np.ndarray,
    #     components: np.ndarray,
    #     tangled_storage: bool = True,
    # ) -> np.ndarray:

    #     components = np.asarray(components, dtype=np.uint32)
    #     self._check_field_components(components)

    #     pass

    def get_field_state_components(
        self,
        components: np.ndarray,
        tangled_storage: bool = True,
    ) -> np.ndarray:
        """Return components of system field state values"""

        components = np.asarray(components, dtype=np.uint32)
        self._check_field_components(components)
        num_components = np.uint32(components.size)

        if num_components == 0:
            return np.array([])

        if num_components == 1:
            tangled_storage = False

        component_field_states = np.zeros((self.num_nodes, num_components))

        # Tangled storage is a switch indicating C-contiguous or
        # F-contigous arrays

        self._cpp_object.get_field_state_components(
            component_field_states,
            num_components,
            components,
            tangled_storage,
        )
        if num_components > 1:
            return component_field_states[:, components]
        return component_field_states
