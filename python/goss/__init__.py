# flake8: noqa
from __future__ import annotations

from . import cli
from . import codegeneration
from . import solvers
from .cli import __author__
from .cli import __license__
from .cli import __version__
from .compilemodule import jit
from .ode import ODE
from .ode import ParameterizedODE
from .systemsolver import ODESystemSolver


def has_openmp() -> bool:
    from . import _gosscpp  # type: ignore

    return _gosscpp.has_openmp()


__all__ = [
    "codegeneration",
    "jit",
    "ODE",
    "solvers",
    "ODESystemSolver",
    "ParameterizedODE",
    "cli",
    "__version__",
]

# If dolfin is present import it
try:
    from . import dolfinutils  # noqa: F401
    from .dolfinutils import DOLFINODESystemSolver  # noqa: F401

    __all__.extend(["dolfinutils", "DOLFINODESystemSolver"])


except ImportError:
    pass
