from typing import NamedTuple
from typing import Optional
import os
import numpy as np
import gotran

from .compilemodule import jit


class LinearizedEval(NamedTuple):
    rhs: np.ndarray
    linear: Optional[np.ndarray] = None


class ODE:
    def __init__(self, *args, **kwargs):

        if isinstance(args[0], gotran.ODE):
            self._cpp_object = jit(args[0])
        elif isinstance(args[0], os.PathLike):
            # Assume this is a gotran ode file
            self._cpp_object = jit(gotran.load_ode(args[0]))

        else:
            # Need to import this module here in order to
            # not cause trouble for cppyy
            from . import _gosscpp

            if isinstance(args[0], _gosscpp.ODE):
                self._cpp_object = args[0]
            else:
                raise ValueError(f"Unknown argument of type {type(args[0])}")

    def copy(self):
        return ODE(self._cpp_object)

    @property
    def num_states(self) -> int:
        """Number of states"""
        return self._cpp_object.num_states()

    @property
    def is_dae(self) -> bool:
        """Return True if the ODE is a
        Differential-algebraic system of equations"""
        return self._cpp_object.is_dae()

    def get_ic(self):
        """Get default initial conditions"""
        return self._cpp_object.get_ic()

    def _check_state_input(self, states) -> np.ndarray:
        msg = (
            "State array has to be of same size as the number"
            f"of states, got {len(states)}, expected {self.num_states}"
        )
        if len(states) != self.num_states:
            raise ValueError(msg)
        if not isinstance(states, np.ndarray):
            return np.array(states)
        return states

    def eval(self, states: np.ndarray, time: float) -> np.ndarray:
        """Evaluate rhs of the ODE"""
        states = self._check_state_input(states)
        values = np.zeros_like(states)
        self._cpp_object.eval(states, time, values)
        return values

    def eval_component(self, component: int, states: np.ndarray, time: float) -> float:
        """Evaluate component idx of the rhs of the ODE"""
        states = self._check_state_input(states)
        return self._cpp_object.eval(component, states, time)

    def compute_jacobian(self, states: np.ndarray, time: float) -> np.ndarray:
        states = self._check_state_input(states)
        jac = np.zeros((states.size, states.size))
        self._cpp_object.compute_jacobian(states, time, jac)
        return jac

    def linearized_eval(
        self, states: np.ndarray, time: float, only_linear: bool = False
    ) -> LinearizedEval:
        states = self._check_state_input(states)
        linear = np.zeros_like(states)
        rhs = np.zeros_like(states)

        self._cpp_object.linearized_eval(states, time, linear, rhs, only_linear)

        if only_linear:
            return LinearizedEval(rhs=rhs)
        return LinearizedEval(linear=linear, rhs=rhs)
