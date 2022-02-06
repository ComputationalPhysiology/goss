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

__all__ = ["jit"]

# System imports
import sys
import dijitso
import hashlib
import pkgconfig
import itertools
from pathlib import Path

# Import Gotran
import gotran
from gotran.common import check_arg, push_log_level, pop_log_level, info, INFO
from gotran.codegeneration.codegenerators import PythonCodeGenerator
from gotran.model.ode import ODE
from gotran.codegeneration.compilemodule import load_module, save_module

# Local imports
from .codegeneration import GossCodeGenerator
from . import cpp


if pkgconfig.exists("goss"):
    goss_pc = pkgconfig.parse("goss")
else:
    raise RuntimeError("Could not find GOSS pkg-config file. Please make sure appropriate paths are set.")


def _jit_generate(class_data, module_name, signature, parameters):
 
    code_c = class_data["code"]
    code_h = ""
    depends = []

    return code_h, code_c, depends

def find_shared_libs(goss_pc):
    shared_libs = []
    suffixes = [".dylib", ".so", ".dll"]

    # Add current directory so that we can support
    # editable installs and virtual environments
    extra_libs = [Path(__file__).absolute().parent, sys.prefix + "/lib"]
    libdirs = goss_pc["library_dirs"] + extra_libs
    for lib in goss_pc["libraries"]:
        for (libdir, suffix) in itertools.product(libdirs, suffixes):
            path = Path(libdir).joinpath(f"lib{lib}").with_suffix(suffix)
            print(path)
            if path.is_file():
                shared_libs.append(path.as_posix())
                break
        else:
            raise RuntimeError(f"Could not find shared library {lib}. Please update your path")
    return shared_libs
    


def jit(ode, field_states=None, field_parameters=None,
        monitored=None, code_params=None, cppargs=None):
    """
    Generate a goss::ODEParameterized from a gotran ode and JIT compile it

    Arguments:
    ----------
    ode : gotran.ODE
        The gotran ode, either as an ODE or as an ODERepresentation
    field_states : list
        A list of state names, which should be treated as field states
    field_parameters : list
        A list of parameter names, which should be treated as field parameters
    monitored : list
        A list of names of intermediates of the ODE. Code for monitoring
        the intermediates will be generated.
    code_params : dict
        Parameters controling the code generation
    cppargs : str
        Default C++ argument passed to the C++ compiler
    """

    # Code generators
    cgen = GossCodeGenerator(ode, field_states, field_parameters, monitored, code_params)
    cgen.params.class_code = True
    pgen = PythonCodeGenerator(cgen.params)
    
    # Create unique module name for this application run
    module_name = "goss_compiled_module_{0}_{1}".format(\
        ode.name, hashlib.sha1((ode.signature() + \
                               repr(code_params) + \
                               repr(field_states) + \
                               repr(field_parameters) + \
                               repr(monitored) + \
                               dijitso.__version__ + \
                               gotran.__version__ + \
                               str(cppargs)).encode()).hexdigest())
    
    # Check cache
    compiled_module = load_module(module_name)

    if compiled_module:
        return getattr(compiled_module, cgen.name)()

    push_log_level(INFO)

    # Init state code
    python_code = pgen.init_states_code(ode)
    cpp_code = cgen.file_code()

    dijitso_params = dijitso.validate_params(dijitso.params.default_params())
    dijitso_params['build']['include_dirs'] = goss_pc["include_dirs"]
    dijitso_params['build']['libs'] = goss_pc["libraries"]
    dijitso_params['build']['lib_dirs'] = goss_pc["library_dirs"]
    # I guess normally dijitso should make this work, but I only managed to get it
    # to work if I used shared libraries. Something to refactor in th future.
    dijitso_params['build']["cxxflags"] = list(dijitso_params['build']["cxxflags"]) + find_shared_libs(goss_pc)
    cpp_data = {"code": cpp_code}

    module, signature = dijitso.jit(
        cpp_data,
        module_name,
        dijitso_params,
        generate=_jit_generate,
    )

    info("Calling GOSS just-in-time (JIT) compiler, this may take some "\
         "time...")
    sys.stdout.flush()

    declaration_form = dict(\
        ModelName = cgen.name, 
        python_code = python_code,
        )

    info(" done")
    pop_log_level()
    sys.stdout.flush()

    # Return an instantiated class
    return getattr(compiled_module, cgen.name)()

