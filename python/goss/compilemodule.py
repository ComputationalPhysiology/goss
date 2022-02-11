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
from pathlib import Path

# Import Gotran
from modelparameters.logger import push_log_level, INFO

# Local imports
from .codegeneration import GossCodeGenerator


here = Path(__file__).parent.absolute()


def cppyy_jit(
    ode,
    field_states=None,
    field_parameters=None,
    monitored=None,
    code_params=None,
):
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
    """

    # Code generators
    cgen = GossCodeGenerator(
        ode,
        field_states,
        field_parameters,
        monitored,
        code_params,
        add_signature_to_name=True,
    )
    cgen.params.class_code = True

    push_log_level(INFO)

    # Init state code
    cpp_code = cgen.file_code()

    import cppyy

    cppyy.add_include_path(here.joinpath("include").as_posix())
    cppyy.add_library_path(here.joinpath("lib").as_posix())
    cppyy.load_library("goss")

    cppyy.cppdef(cpp_code)

    _cppyygbl = __import__("cppyy.gbl", fromlist=[f"create_{cgen.name}"])
    submodule = getattr(_cppyygbl, f"create_{cgen.name}")()

    import cppyy.ll

    return cppyy.ll.as_ctypes(submodule)


def make_ode(python_object):
    from . import _gosscpp

    return _gosscpp.make_ode(python_object.value)


def jit(
    ode,
    field_states=None,
    field_parameters=None,
    monitored=None,
    code_params=None,
):
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
    """

    python_object = cppyy_jit(
        ode=ode,
        field_states=field_states,
        field_parameters=field_parameters,
        monitored=monitored,
        code_params=code_params,
    )
    return make_ode(python_object=python_object)
