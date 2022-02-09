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

    def copy(self):
        return self.__class__(self._cpp_object)

    def forward(self, y: np.ndarray, t: float, interval: float) -> np.ndarray:
        self._cpp_object.forward(y, t, interval)
        return y

    @property
    def num_states(self):
        return self._cpp_object.num_states()


class ExplicitEuler(ODESolver):
    pass


class RL1(ODESolver):
    pass


class ImplicitODESolver(ODESolver):
    pass
