"""Expose names from the compiled cpp modules to the goss.cpp namespace"""
# Copyright (C) 2012 Johan Hake
#
# This file is part of GOSS.
#
# GOSS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GOSS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GOSS. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2012-10-11
# Last changed: 2012-10-11
# Import pure python modules
from __future__ import annotations

from . import codegeneration
from . import solvers
from .compilemodule import jit
from .ode import ODE
from .systemsolver import ODESystemSolver

goss_solvers = [solvers.RL1]

__all__ = ["codegeneration", "jit", "ODE", "goss_solvers", "ODESystemSolver"]

# If dolfin is present import it
# try:
#     from . import dolfinutils

#     __all__.extend(dolfinutils.__all__)
# except ModuleNotFoundError as e:
#     pass
