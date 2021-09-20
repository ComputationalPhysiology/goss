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
import os
import re
import numpy
import instant
import hashlib
import types

# Import Gotran
import gotran
from gotran.common import check_arg, push_log_level, pop_log_level, info, INFO
from gotran.codegeneration.codegenerators import PythonCodeGenerator
from gotran.model.ode import ODE

# Local imports
from .codegeneration import GossCodeGenerator
from . import cpp

# Set log level of instant
instant.set_log_level("WARNING")

_additional_declarations = r"""
%init%{{
import_array();
%}}

%include <exception.i>
%feature("autodoc", "1");
%include <std_string.i>
%include <goss/swig/typemaps.i>
%include <goss/swig/exceptions.i>

%include <boost_shared_ptr.i>

%shared_ptr(goss::ODE)
%shared_ptr(goss::ParameterizedODE)
%shared_ptr(goss::{ModelName})

// GOSS Import statements
%import(module="goss.cpp") "goss/ODE.h"
%import(module="goss.cpp") "goss/ParameterizedODE.h"

// Modifications of the wrapped C++ class
%ignore goss::{ModelName}::get_ic;
%extend goss::{ModelName}{{
%pythoncode%{{
{python_code}
%}}
}}

// Rename 
%rename(_eval) goss::{ModelName}::eval(const double* states, double time, double* values);
%rename(eval_component) goss::{ModelName}::eval(uint, const double*, double);
"""

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
                               #instant.get_swig_version() + \
                               instant.__version__ + \
                               gotran.__version__ + \
                               str(cppargs)).encode()).hexdigest())
    
    # Check cache
    compiled_module = instant.import_module(module_name)

    if compiled_module:
        return getattr(compiled_module, cgen.name)()

    push_log_level(INFO)

    # Init state code
    python_code = pgen.init_states_code(ode)
    cpp_code = cgen.class_code()

    info("Calling GOSS just-in-time (JIT) compiler, this may take some "\
         "time...")
    sys.stdout.flush()

    # Configure instant and add additional system headers
    instant_kwargs = configure_instant()

    instant_kwargs["cppargs"] = cppargs or instant_kwargs["cppargs"]
    instant_kwargs["cmake_packages"] = ["GOSS"] 

    declaration_form = dict(\
        ModelName = cgen.name, 
        python_code = python_code,
        )

    # Compile extension module with instant
    compiled_module = instant.build_module(\
        code = cpp_code,
        additional_declarations = _additional_declarations.format(\
            **declaration_form),
        signature = module_name,
        **instant_kwargs)

    info(" done")
    pop_log_level()
    sys.stdout.flush()

    # Return an instantiated class
    return getattr(compiled_module, cgen.name)()

def configure_instant():
    """
    Check system requirements

    Returns a dict with kwargs that can be passed to instant.build_module.
    """
    instant_kwargs = {}
    swig_include_dirs = []

    # Let swig see the installed gillstep swig files
    swig_include_dirs = []
    goss_include_found = False

    # Check that the form compiler will use the same swig version
    # that PyGOSS was compiled with
    if not instant.check_swig_version(cpp.__swigversion__, same=True):
        raise OSError("""GOSS was not compiled with the present version of swig.
Install swig version {0} or recompiled GOSS with present swig
""".format(cpp.__swigversion__))

    instant_kwargs['system_headers'] = ["boost/shared_ptr.hpp",
                                        "boost/make_shared.hpp",
                                        "cmath",
                                        "stdexcept",
                                        "numpy/arrayobject.h",
                                        "goss/ParameterizedODE.h",
                                        "goss/Timer.h",
                                        ]
    instant_kwargs['swigargs'] =['-O -c++']
    instant_kwargs['cppargs'] = ['-O2']

    return instant_kwargs
