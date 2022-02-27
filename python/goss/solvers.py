from __future__ import annotations

import abc
from enum import Enum
from typing import Any
from typing import Optional

import numpy as np

from .ode import ODE


class ODESolver(abc.ABC):
    """Base class for ODESolver.

    Your can instantiate an ODE solver in many different ways,
    but the general theme is that an ode solver need an ode-model
    that is typically coming from a gotran ode-file.

    The default way to instantiate a solver is to provide the
    ode as the first argument

    .. code:: python

        solver = ExplicitEuler(ode)

    where `ode` is of type `goss.ode.ODE`.
    Another way is to create an empty ODESolver

    .. code:: python

        solver = ExplicitEuler()

    but then you need to attach an ode later in order to use it

    .. code:: python

        solver.attach(ode)

    """

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
        elif isinstance(args[0], _gosscpp.ODESolver):
            self._cpp_object = args[0]
        else:
            raise ValueError(f"Unknown argument of type {type(args[0])}")

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"

    def reset(self) -> None:
        self._cpp_object.reset()

    def attach(self, ode: ODE) -> None:
        self._cpp_object.attach(ode._cpp_object)

    @property
    def is_adaptive(self) -> bool:
        """Flag to indicate whether the solver is adaptive for not"""
        return self._cpp_object.is_adaptive()

    @staticmethod
    def default_parameters() -> dict[str, Any]:
        """Dictionary with the default parameters"""
        return {}

    def update_parameters(self, parameters: dict[str, Any]):
        """Update the parameters given a dictionary

        Parameters
        ----------
        parameters : dict[str, Any]
            The new parameters
        """
        for k, v in parameters.items():
            self.set_parameter(k, v)

    @property
    def internal_time_step(self) -> float:
        """Time step used internally by the solver"""
        return self._cpp_object.get_internal_time_step()

    @internal_time_step.setter
    def internal_time_step(self, time_step: float) -> None:
        self._cpp_object.set_internal_time_step(time_step)

    @property
    def parameters(self) -> dict[str, Any]:
        """Dictionaty with the current parameters"""
        parameters = self.__class__.default_parameters()
        return {name: self.get_parameter(name) for name in parameters}

    def get_parameter(self, name: str) -> Any:
        """Get the current value of a parameter

        Parameters
        ----------
        name : str
            The name of the paramrter of interest

        Returns
        -------
        Any
            The value of the parameter
        """
        return getattr(self._cpp_object, name)

    def set_parameter(self, name: str, value: Any):
        """Set the value of a parameter

        Parameters
        ----------
        name : str
            Name of the parameter you want to set
        value : Any
            The new value of the parameter

        Raises
        ------
        KeyError
            If the name is not a valid parameter
        TypeError
            If the type of the new value does not match the original value
        """
        parameters = self.__class__.default_parameters()
        if name not in parameters:
            raise KeyError(
                f"Invalid parameter {name}, expected one of {parameters.keys()}",
            )
        if not isinstance(value, type(parameters[name])):
            raise TypeError(
                (
                    f"Expected parameter {name} to be of type"
                    f"{type(parameters[name])}, got {type(value)}"
                ),
            )
        setattr(self._cpp_object, name, value)

    def get_ode(self) -> ODE:
        """Get the ODE from the cpp object"""
        return ODE(self._cpp_object.get_ode())

    def copy(self):
        """Make a copy of the solver"""
        return self.__class__(self._cpp_object.copy())

    def forward(self, y: np.ndarray, t: float, interval: float):
        """Do one forward iteration

        Parameters
        ----------
        y : np.ndarray
            Array with the current state values. Note that
            this array will be mutatated and contain the updated
            states after call to this function.
        t : float
            Time point
        interval : float
            Interval to step
        """
        self._cpp_object.forward(y, t, interval)

    def solve(
        self,
        t: np.ndarray,
        y0: Optional[np.ndarray] = None,
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
        if y0 is None:
            # Use the initial conditions from the ODE
            y0 = self.get_ode().get_ic()
        num_steps = t.size
        y = np.zeros((num_steps, self.num_states))

        y[0, :] = y0
        self._cpp_object.solve(y, y0, t, num_steps)
        return y

    @property
    def num_states(self):
        """Number of states in the ODE"""
        return self._cpp_object.num_states()


class ExplicitEuler(ODESolver):
    pass


class RL1(ODESolver):
    pass


class RL2(ODESolver):
    pass


class RK2(ODESolver):
    pass


class RK4(ODESolver):
    pass


class GRL1(ODESolver):
    @staticmethod
    def default_parameters() -> dict[str, Any]:
        names = ODESolver.default_parameters()
        names.update({"delta": 1e-8})
        return names


class GRL2(GRL1):
    pass


class ImplicitODESolver(ODESolver, abc.ABC):
    @staticmethod
    def default_parameters() -> dict[str, Any]:
        names = ODESolver.default_parameters()
        names.update(
            {
                "eta_0": 1.0,
                "kappa": 0.1,
                "relative_tolerance": 1e-12,
                "max_iterations": 30,
                "max_relative_previous_residual": 0.01,
                "always_recompute_jacobian": False,
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
        dt: float,
        alpha: float,
    ) -> None:
        return self._cpp_object.compute_factorized_jacobian(y, t, dt, alpha)


class BasicImplicitEuler(ImplicitODESolver):
    """FIXME: This is not working as expected"""

    pass


class ImplicitEuler(ImplicitODESolver):
    @staticmethod
    def default_parameters() -> dict[str, Any]:
        names = ImplicitODESolver.default_parameters()
        names.update(
            {
                "num_refinements_without_always_recomputing_jacobian": 2,
                "min_dt": 0.0001,
            },
        )
        return names


class ThetaSolver(ImplicitODESolver):
    @staticmethod
    def default_parameters() -> dict[str, Any]:
        names = ImplicitODESolver.default_parameters()
        names.update(
            {
                "num_refinements_without_always_recomputing_jacobian": 2,
                "min_dt": 0.0001,
                "theta": 0.5,
            },
        )
        return names


class AdaptiveSolver(ODESolver, abc.ABC):
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


class AdaptiveExplicitSolver(AdaptiveSolver, abc.ABC):
    pass


class RKF32(AdaptiveExplicitSolver):
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


class AdaptiveImplicitSolver(ImplicitODESolver, AdaptiveSolver, abc.ABC):
    pass


class ESDIRK(AdaptiveImplicitSolver):
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


class ESDIRK4O32(ESDIRK):
    """FIXME: This is not working as expected"""

    pass


class ESDIRK23a(ESDIRK):
    @staticmethod
    def default_parameters() -> dict[str, Any]:
        names = AdaptiveImplicitSolver.default_parameters()
        names.update(
            {
                "num_refinements_without_always_recomputing_jacobian": 2,
                "min_dt": 0.001,
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


class GOSSSolvers(Enum):
    ExplicitEuler = "ExplicitEuler"
    RK2 = "RK2"
    RK4 = "RK4"
    RL1 = "RL1"
    RL2 = "RL2"
    GRL1 = "GRL1"
    GRL2 = "GRL2"
    # BasicImplicitEuler = "BasicImplicitEuler"
    ImplicitEuler = "ImplicitEuler"
    ThetaSolver = "ThetaSolver"
    RKF32 = "RKF32"
    ESDIRK23a = "ESDIRK23a"
    # ESDIRK4O32 = "ESDIRK4O32"


class GOSSImplicitSolvers(Enum):
    # BasicImplicitEuler = "BasicImplicitEuler"
    ImplicitEuler = "ImplicitEuler"
    ThetaSolver = "ThetaSolver"
    ESDIRK23a = "ESDIRK23a"
    # ESDIRK4O32 = "ESDIRK4O32"


class GOSSExplicitSolvers(Enum):
    ExplicitEuler = "ExplicitEuler"
    RK2 = "RK2"
    RK4 = "RK4"
    RL1 = "RL1"
    RL2 = "RL2"
    GRL1 = "GRL1"
    GRL2 = "GRL2"
    RKF32 = "RKF32"


class GOSSNonAdaptiveSolvers(Enum):
    ExplicitEuler = "ExplicitEuler"
    RK2 = "RK2"
    RK4 = "RK4"
    RL1 = "RL1"
    RL2 = "RL2"
    GRL1 = "GRL1"
    GRL2 = "GRL2"
    # BasicImplicitEuler = "BasicImplicitEuler"
    ImplicitEuler = "ImplicitEuler"
    ThetaSolver = "ThetaSolver"


class GOSSAdaptiveSolvers(Enum):
    RKF32 = "RKF32"
    ESDIRK23a = "ESDIRK23a"
    # ESDIRK4O32 = "ESDIRK4O32"


solver_mapper = {
    "ExplicitEuler": ExplicitEuler,
    "RK2": RK2,
    "RK4": RK4,
    "RL1": RL1,
    "RL2": RL2,
    "GRL1": GRL1,
    "GRL2": GRL2,
    # "BasicImplicitEuler": BasicImplicitEuler,
    "ImplicitEuler": ImplicitEuler,
    "ThetaSolver": ThetaSolver,
    "RKF32": RKF32,
    "ESDIRK23a": ESDIRK23a,
    # "ESDIRK4O32": ESDIRK4O32,
}
