from __future__ import annotations

from typing import Any

import numpy as np

from .ode import ODE


class ODESolver:
    def __init__(self, *args, **kwargs):
        from . import _gosscpp

        # Get the name of the class
        name = self.__class__.__name__
        # TODO: Add a check to make sure that the name is valid

        # Grab the C++ constructor
        cpp_solver = getattr(_gosscpp, name)

        # TODO: Implement a protocol for the cpp_object
        if len(args) == 0:
            self._cpp_object = cpp_solver()
        elif isinstance(args[0], ODE):
            self._cpp_object = cpp_solver(args[0]._cpp_object)
        else:
            raise ValueError(f"Unknown argument of type {type(args[0])}")

    @property
    def is_adaptive(self) -> bool:
        return self._cpp_object.is_adaptive()

    @property
    def parameter_names(self) -> list[str]:
        return ["ldt"]

    @property
    def parameters(self):
        return {name: getattr(self._cpp_object, name) for name in self.parameter_names}

    def set_parameter(self, name: str, value: Any):
        assert name in self.parameter_names
        setattr(self._cpp_object, name, value)

    def get_internal_time_step(self) -> float:
        # Don'r really know what this is?
        return self._cpp_object.get_internal_time_step()

    def get_ode(self) -> ODE:
        return ODE(self._cpp_object.get_ode())

    def copy(self):
        return self.__class__(self._cpp_object)

    def forward(self, y: np.ndarray, t: float, interval: float):
        # FIXME: Consider making this pure
        self._cpp_object.forward(y, t, interval)

    @property
    def num_states(self):
        return self._cpp_object.num_states()


class ExplicitEuler(ODESolver):
    pass


class RL1(ODESolver):
    pass


class ImplicitODESolver(ODESolver):
    pass
