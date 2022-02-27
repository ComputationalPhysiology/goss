from __future__ import annotations

import ctypes
import typing
from pathlib import Path

import gotran
from modelparameters.logger import INFO
from modelparameters.logger import push_log_level

from .codegeneration import GossCodeGenerator

# Import Gotran
# Local imports


here = Path(__file__).parent.absolute()


def cppyy_jit(
    ode: gotran.ODE,
    field_states: list[str] = None,
    field_parameters: list[str] = None,
    monitored: list[str] = None,
    code_params: dict[str, typing.Any] = None,
) -> ctypes.c_void_p:
    """
    Generate a goss::ODEParameterized from a gotran ode and JIT compile it

    Parameters
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
    return code_to_submodule(cpp_code, cgen.name)


def code_to_submodule(code: str, name):

    import cppyy

    cppyy.add_include_path(here.joinpath("include").as_posix())
    cppyy.add_library_path(here.joinpath("lib").as_posix())
    cppyy.load_library("goss")

    cppyy.cppdef(code)

    _cppyygbl = __import__("cppyy.gbl", fromlist=[f"create_{name}"])
    submodule = getattr(_cppyygbl, f"create_{name}")()

    import cppyy.ll

    return cppyy.ll.as_ctypes(submodule)


def make_ode(python_object):
    from . import _gosscpp

    return _gosscpp.make_ode(python_object.value)


def jit(
    ode: gotran.ODE,
    field_states: list[str] = None,
    field_parameters: list[str] = None,
    monitored: list[str] = None,
    code_params: dict[str, typing.Any] = None,
):
    """
    Generate a goss::ODEParameterized from a gotran ode and JIT compile it

    Parameters
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
