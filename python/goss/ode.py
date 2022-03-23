from __future__ import annotations

import os
from pathlib import Path
from typing import NamedTuple
from typing import Optional

import gotran
import numpy as np

from . import compilemodule


class LinearizedEval(NamedTuple):
    rhs: np.ndarray
    linear: Optional[np.ndarray] = None


def load_file(path: Path, **kwargs):

    while not path.is_file():
        if path.suffix == "":
            # Assume extension is ".ode"
            path = path.with_suffix(".ode")
        else:
            raise FileNotFoundError(f"File {path} does not exist")

    if path.suffix == ".ode":
        return compilemodule.jit(gotran.load_ode(path), **kwargs)
    elif path.suffix in [".h", ".hpp"]:
        name = kwargs.get("name", "")
        if name == "":
            raise ValueError("Missing name of cpp module")
        with open(path, "r") as f:
            code = f.read()
        return compilemodule.make_ode(compilemodule.code_to_submodule(code, name))

    raise RuntimeError(f"Invalid extension {path.suffix}")


class ODE:
    def __init__(self, *args, **kwargs):

        if isinstance(args[0], gotran.ODE):
            self._cpp_object = compilemodule.jit(args[0], **kwargs)
        elif isinstance(args[0], (os.PathLike, str)):
            # Assume this is a gotran ode file
            path = Path(args[0])
            self._cpp_object = load_file(path, **kwargs)

        else:
            # Need to import this module here in order to
            # not cause trouble for cppyy
            from . import _gosscpp

            if isinstance(args[0], (_gosscpp.ODE, _gosscpp.ParameterizedODE)):
                self._cpp_object = args[0]
            else:
                raise ValueError(f"Unknown argument of type {type(args[0])}")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._cpp_object})"

    def __str__(self) -> str:
        return f"{self.__class__.__name__} with {self.num_states} states"

    def copy(self):
        return self.__class__(self._cpp_object)

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
        """Compute numerical jacobian"""
        states = self._check_state_input(states)
        jac = np.zeros((states.size, states.size))
        self._cpp_object.compute_jacobian(states, time, jac)
        return jac

    def linearized_eval(
        self,
        states: np.ndarray,
        time: float,
        only_linear: bool = False,
    ) -> LinearizedEval:
        """Evaluate the linearized rhs"""
        states = self._check_state_input(states)
        linear = np.zeros_like(states)
        rhs = np.zeros_like(states)

        self._cpp_object.linearized_eval(states, time, linear, rhs, only_linear)

        if only_linear:
            return LinearizedEval(rhs=rhs)
        return LinearizedEval(linear=linear, rhs=rhs)

    def lu_factorize(self, mat: np.ndarray) -> None:
        """In place LU Factorize matrix (jacobian)"""
        self._cpp_object.lu_factorize(mat)

    def forward_backward_subst(
        self,
        mat: np.ndarray,
        b: np.ndarray,
        x: np.ndarray,
    ) -> None:
        """Forward/Backward substitution of factorized matrix"""
        self._cpp_object.forward_backward_subst(mat, b, x)


class ParameterizedODE(ODE):
    def __init__(self, *args, **kwargs):
        """FIXME: Write about all possible ways
        to initialize object
        """
        if "num_states" in kwargs and "num_parameters" in kwargs:
            num_states = kwargs.get("num_states")
            num_parameters = kwargs.get("num_parameters")
            num_field_states = kwargs.get("num_field_states", 0)
            num_field_parameters = kwargs.get("num_field_parameters", 0)
            num_monitored = kwargs.get("num_monitored", 0)
            from . import _gosscpp

            self._cpp_object = _gosscpp.ParameterizedODE(
                num_states,
                num_parameters,
                num_field_states,
                num_field_parameters,
                num_monitored,
            )
        else:
            super().__init__(*args, **kwargs)

    @property
    def num_parameters(self) -> int:
        return self._cpp_object.num_parameters()

    @property
    def num_field_parameters(self) -> int:
        return self._cpp_object.num_field_parameters()

    @property
    def num_field_states(self) -> int:
        return self._cpp_object.num_field_states()

    @property
    def num_monitored(self) -> int:
        return self._cpp_object.num_monitored()

    def eval_monitored(self, states: np.ndarray, time: float) -> np.ndarray:
        assert len(states) == self.num_states
        monitored = np.zeros(self.num_monitored)
        self._cpp_object.eval_monitored(states, time, monitored)
        return monitored

    def monitored_values(self, states: np.ndarray, time: np.ndarray) -> np.ndarray:
        num_time_steps = len(time)
        assert states.shape == (num_time_steps, self.num_states)
        monitored = np.zeros((num_time_steps, self.num_monitored))
        m = np.zeros(self.num_monitored)
        self._cpp_object.monitored_values(states, time, monitored, m)
        return monitored

    def set_parameter(self, name: str, value: float) -> None:
        self._cpp_object.set_parameter(name, value)

    def get_parameter(self, name: str) -> float:
        return self._cpp_object.get_parameter(name)

    @property
    def parameters(self) -> dict[str, float]:
        # FIXME: This can potentially be slow if we have a lot of parameters
        return {name: self.get_parameter(name) for name in self.parameter_names}

    @property
    def state_names(self) -> list[str]:
        return self._cpp_object.get_state_names()

    @property
    def field_state_names(self) -> list[str]:
        return self._cpp_object.get_field_state_names()

    @property
    def parameter_names(self) -> list[str]:
        return self._cpp_object.get_parameter_names()

    @property
    def field_parameter_names(self) -> list[str]:
        return self._cpp_object.get_field_parameter_names()

    @property
    def field_state_indices(self) -> list[str]:
        return self._cpp_object.get_field_state_indices()

    @property
    def monitored_names(self) -> list[str]:
        return self._cpp_object.get_monitored_names()
