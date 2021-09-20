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

__all__ = []

# Import pure python modules
from . import codegeneration as codegeneration
from . import compilemodule as compilemodule

from .codegeneration import *
from .compilemodule import *

__all__.extend(codegeneration.__all__)
__all__.extend(compilemodule.__all__)

#--- Imports the SWIG-generated Python code (C++ interface) ---

_not_includes = ["new_instancemethod", "SHARED_PTR_DISOWN", \
                 "weakref", "weakref_proxy", "cvar"]

from . import cpp

# Add the module to the global namespace
__all__.append("cpp")

# Collect all classes to outrule static class methods
_all_classes = []

# Iterate over the modules attributes and add names to global namespace
for name in sorted(cpp.__dict__.keys()):

    # Skip some includes
    if name in _not_includes:
        continue

    # Do not add an attrbute name which contains "swigregister" or
    # starts with "_"
    if "swigregister" in name or name[0] == "_":
        continue

    # Get attr
    attr = getattr(cpp, name)

    # Check for class
    if isinstance(attr, type):
        _all_classes.append(name)
    else:
        # If not a class, check if first part of function name is the
        # same as a class name
        if name.split("_")[0] in _all_classes:
            continue
        
    # Add the attribute to the global namespace
    globals()[name] = attr
    __all__.append(name)

# Add version 
# -----------
globals()["__version__"] = getattr(cpp, "__gossversion__")
globals()["__swigversion__"] = getattr(cpp, "__swigversion__")

# Get names of all valid goss.ODESolvers
abstract_solvers = ["ImplicitODESolver", "AdaptiveImplicitSolver", "AdaptiveExplicitSolver"]
goss_solvers = [name for name, attr in cpp.__dict__.items() \
                 if isinstance(attr, type) and issubclass(attr, cpp.ODESolver) \
                and name not in abstract_solvers]

goss_solvers.remove("ODESolver")
__all__.append("goss_solvers")

# If dolfin is present import it
try:
    from . import dolfinutils
    from .dolfinutils import *
    __all__.extend(dolfinutils.__all__)
except Exception as e:
    pass
    

# Set debug level
DEBUG = DBG

# Clean up
del name, attr