from __future__ import annotations

import abc
from collections import namedtuple
from typing import Any
from typing import Optional

import numpy as np

from .ode import ODE

Solution = namedtuple("Solution", ["y", "t"])


class ODESolver(abc.ABC):
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
    def parameter_names(self) -> dict[str, Any]:
        return {"ldt": float}

    @property
    def parameters(self):
        return {name: getattr(self._cpp_object, name) for name in self.parameter_names}

    def set_parameter(self, name: str, value: Any):
        assert name in self.parameter_names
        if not isinstance(value, self.parameter_names[name]):
            raise TypeError(
                (
                    f"Expected parameter {name} to be of type"
                    f"{self.parameter_names[name]}, got {type(value)}"
                ),
            )
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

    def solve(
        self,
        start: float,
        end: float,
        dt: float,
        y0: Optional[np.ndarray] = None,
        skip_n: int = 1,
    ) -> Solution:
        if y0 is None:
            # Use the initial conditions from the ODE
            y0 = self.get_ode().get_ic()
        interval = end - start
        num_steps = int(np.ceil(interval / dt))
        y = np.zeros((num_steps + 1, self.num_states))
        t = np.arange(start, end + dt, dt)
        y[0, :] = y0
        self._cpp_object.solve(y, y0, t, num_steps, skip_n)
        return Solution(y, t)

    @property
    def num_states(self):
        return self._cpp_object.num_states()


class ExplicitEuler(ODESolver):
    pass


class RL1(ODESolver):
    pass


class GRL1(ODESolver):
    @property
    def parameter_names(self) -> dict[str, Any]:
        names = super().parameter_names
        names.update({"delta": float})
        return names


class ImplicitODESolver(ODESolver, abc.ABC):
    @property
    def parameter_names(self) -> dict[str, Any]:
        names = super().parameter_names
        names.update(
            {
                "kappa": float,
                "relative_tolerance": float,
                "max_iterations": int,
                "max_relative_previous_residual": float,
                "always_recompute_jacobian": bool,
            },
        )
        return names

    def num_jac_comp(self) -> int:
        """Return the number of times the jacobian
        has been computed"""
        return self._cpp_object.num_jac_comp()

    def compute_factorized_jacobian(
        self,
        y: np.ndarray,
        t: float,
        alpha: float,
    ) -> None:
        self._cpp_object.compute_factorized_jacobian(y, t, alpha)


class ThetaSolver(ImplicitODESolver):
    @property
    def parameter_names(self) -> dict[str, Any]:
        names = super().parameter_names
        names.update(
            {
                "num_refinements_without_always_recomputing_jacobian": int,
                "min_dt": float,
                "theta": float,
            },
        )
        return names


class AdaptiveImplicitSolver(ImplicitODESolver, abc.ABC):
    def get_current_time(self) -> float:
        return self._cpp_object.get_current_time()

    def get_current_time_step(self) -> float:
        return self._cpp_object.get_current_time_step()

    def get_num_accepted(self) -> int:
        return self._cpp_object.get_num_accepted()

    def get_num_rejected(self) -> int:
        return self._cpp_object.get_num_rejected()

    def set_single_step_mode(self, mode: bool) -> None:
        self._cpp_object.set_single_step_mode(mode)

    def set_tol(self, atol: float, rtol: float = 1e-8) -> None:
        self._cpp_object.set_tol(atol, rtol)

    def set_iord(self, iord: int) -> None:
        self._cpp_object.set_iord(iord)

    @property
    def atol(self):
        return self._cpp_object.get_atol()

    @property
    def rtol(self):
        return self._cpp_object.get_rtol()

    @property
    def iord(self) -> int:
        return self._cpp_object.get_iord()


class ESDIRK23a(AdaptiveImplicitSolver):
    @property
    def parameter_names(self) -> dict[str, Any]:
        names = super().parameter_names
        names.update(
            {
                "num_refinements_without_always_recomputing_jacobian": int,
                "min_dt": float,
            },
        )
        return names

    @property
    def nfevals(self):
        """Number of right hand side evaluations"""
        return self._cpp_object.nfevals

    @property
    def ndtsa(self):
        """Number of accepted timesteps"""
        return self._cpp_object.ndtsa

    @property
    def ndtsr(self):
        """Number of rejected timesteps"""
        return self._cpp_object.ndtsr
