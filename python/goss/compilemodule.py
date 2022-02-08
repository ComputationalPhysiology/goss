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

# Import Gotran
from modelparameters.logger import push_log_level, INFO
# Local imports
from .codegeneration import GossCodeGenerator
from . import _gosscpp as cpp



def jit(
    ode,
    field_states=None,
    field_parameters=None,
    monitored=None,
    code_params=None,
    cppargs=None,
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
    cppargs : str
        Default C++ argument passed to the C++ compiler
    """

    # Code generators
    cgen = GossCodeGenerator(
        ode, field_states, field_parameters, monitored, code_params
    )
    cgen.params.class_code = True
 
    push_log_level(INFO)

    # Init state code
    cpp_code = cgen.file_code()

    import cppyy.ll

    cppyy.cppdef(cpp_code)
    from cppyy.gbl import create_ODE

    submodule = create_ODE()
    python_object = cppyy.ll.as_ctypes(submodule)

    return cpp.make_ode(python_object.value)
   