from __future__ import annotations

try:
    from importlib.metadata import metadata
except ImportError:
    from importlib_metadata import metadata  # type:ignore

from . import cli
from . import codegeneration
from . import solvers
from .compilemodule import jit
from .ode import ODE
from .ode import ParameterizedODE
from .systemsolver import ODESystemSolver


def has_openmp() -> bool:
    from . import _gosscpp  # type: ignore

    return _gosscpp.has_openmp()


meta = metadata("pygoss")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]


__all__ = [
    "codegeneration",
    "jit",
    "ODE",
    "solvers",
    "ODESystemSolver",
    "ParameterizedODE",
    "cli",
]

# If dolfin is present import it
try:
    from . import dolfinutils  # noqa: F401
    from .dolfinutils import DOLFINODESystemSolver  # noqa: F401

    __all__.extend(["dolfinutils", "DOLFINODESystemSolver"])


except ImportError:
    pass
