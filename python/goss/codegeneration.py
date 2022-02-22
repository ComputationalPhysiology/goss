from __future__ import annotations

import hashlib
from enum import Enum
from typing import Any

import gotran
from gotran.codegeneration.algorithmcomponents import factorized_jacobian_expressions
from gotran.codegeneration.algorithmcomponents import forward_backward_subst_expressions
from gotran.codegeneration.algorithmcomponents import jacobian_expressions
from gotran.codegeneration.algorithmcomponents import linearized_derivatives
from gotran.codegeneration.algorithmcomponents import monitored_expressions
from gotran.codegeneration.algorithmcomponents import rhs_expressions
from gotran.codegeneration.codegenerators import CppCodeGenerator
from gotran.common.options import parameters
from gotran.model.expressions import Expression
from gotran.model.ode import ODE
from modelparameters.logger import error
from modelparameters.utils import check_arg
from modelparameters.utils import check_kwarg
from pydantic import BaseModel
from pydantic import Field

# Gotran imports


_file_template = """#ifndef {MODELNAME}_H_IS_INCLUDED
#define {MODELNAME}_H_IS_INCLUDED
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <memory>
#include <stdexcept>
#include <cmath>

#include <goss/Timer.h>
#include <goss/ParameterizedODE.h>

{CLASS_DECLARATION}
#endif
"""

_class_template = """namespace goss {{

  // Implementation of gotran generated ODE
  class {ModelName} : public ParameterizedODE
  {{
  public:

    // Constructor
    {ModelName}() : ParameterizedODE({num_states}, {num_parameters}, {num_field_states}, {num_field_parameters}, {num_monitored}){variable_initialization}

    {{
{constructor}{constructor_linear_terms}
    }}

    // Copy constructor
    {ModelName}(const {ModelName}& ode) : ParameterizedODE(ode){variable_initialization_copy_constructor}
    {{
      // Do nothing
    }}

    // Evaluate rhs of the ODE
    void eval(const double* states, double time, double* values)
    {{

      //Timer timer_(\"Evaluation of rhs\");
{eval_code}
    }}

{jacobian_code}{factorizing_code}{fb_substitution_code}{linearized_eval_code}{eval_componentwise_code}
    // Get default initial conditions
    void get_ic(goss::DoubleVector *values) const
    {{
{initial_condition_code}
    }}

    // Return a copy of the ODE
    std::shared_ptr<ODE> copy() const
    {{
      return std::make_shared<{ModelName}>(*this);
    }}

    // Evaluate the monitored intermediates
    void eval_monitored(const double* states, double time, double* monitored) const
    {{

      //Timer timer_(\"Evaluation of monitored.\");
{monitored_evaluation_code}
    }}

    // Set all field parameters
    void set_field_parameters(const double* field_params)
    {{
{set_field_parameters_code}
    }}

  private:
{variable_declaration}

  }};

}}

extern "C" DLL_EXPORT goss::ParameterizedODE * create_{ModelName}()
{{
  return new goss::{ModelName};
}}

"""

_no_monitored_snippet = """\n      // No monitored
      throw std::runtime_error(\"No monitored in the \\'{0}\\' model.\");"""

_file_form = dict(
    MODELNAME="NOT_IMPLEMENTED",
    CLASS_DECLARATION="NOT_IMPLEMENTED",
)

_class_form = dict(
    ModelName="NOT_IMPLEMENTED",
    linearized_base_initialization="",
    num_states="NOT_IMPLEMENTED",
    num_parameters=0,
    num_field_states=0,
    num_field_parameters=0,
    num_monitored=0,
    state_names_ctr="NOT_IMPLEMENTED",
    variable_initialization="NOT_IMPLEMENTED",
    variable_initialization_copy_constructor="NOT_IMPLEMENTED",
    constructor="",
    constructor_linear_terms="",
    eval_code="NOT_IMPLEMENTED",
    initial_condition_code="NOT_IMPLEMENTED",
    monitored_evaluation_code="",
    set_field_parameters_code="",
    variable_declaration="NOT_IMPLEMENTED",
    eval_componentwise_code="",
    linearized_eval_code="",
    jacobian_code="",
    factorizing_code="",
    fb_substitution_code="",
)


class StateRepr(Enum):
    named = "named"
    array = "array"


class BodyRepr(Enum):
    named = "named"
    array = "array"
    reused_array = "reused_array"


class OptimizeExprs(Enum):
    none = "none"
    numerals = "numerals"
    numerals_symbols = "numerals_symbols"


class GossCodeGeneratorParameters(BaseModel):
    state_repr: StateRepr = Field(
        StateRepr.named,
        description="Representation of the state",
        type="|".join(StateRepr._member_names_),
    )
    body_repr: BodyRepr = Field(
        BodyRepr.named,
        description="Representation of the body",
        type="|".join(BodyRepr._member_names_),
    )
    use_cse: bool = Field(
        False,
        description="Use sympy cse to optimize common sub expressions",
    )
    generate_forward_backward_subst: bool = Field(
        False,
        description="Generate forward backward substitions for a factorized jacobian",
    )
    generate_jacobian: bool = Field(False, description="Generate Jacobian matrix")
    generate_lu_factorization: bool = Field(
        False,
        description="Geneate Lu factorization",
    )
    optimize_exprs: OptimizeExprs = Field(
        OptimizeExprs.none,
        description="Optimize expressions",
        type="|".join(OptimizeExprs._member_names_),
    )

    def update(self, params):
        # Hack to make this mimic the ParameterDict
        params = params or {}
        self.__init__(**params)

    @classmethod
    def print_defaults(cls):
        from rich.console import Console
        from rich.table import Table

        schema = cls.schema(ref_template="{model}")

        table = Table(title=schema["title"])
        table.add_column("Parameters", justify="right", style="cyan", no_wrap=True)
        table.add_column("Default value", style="magenta")
        table.add_column("Type", style="yellow")
        table.add_column("Description", style="green")

        for k, v in schema["properties"].items():

            t = v.get("type", "")
            desc = v.get("description", "")

            table.add_row(k, repr(v["default"]), t, desc)
        console = Console()
        console.print(table)


class GossCodeGenerator(CppCodeGenerator):
    """
    Class for generating an implementation of a goss ODE
    """

    @staticmethod
    def default_parameters():
        return GossCodeGeneratorParameters()

    def __init__(
        self,
        ode,
        field_states: list[str] = None,
        field_parameters: list[str] = None,
        monitored: list[str] = None,
        code_params: dict[str, Any] = None,
        add_signature_to_name: bool = False,
    ):
        """
        A class for generating a C++ subclass of a goss::ODEParameterized

        Arguments:
        ----------
        ode : gotran.ODE
            The gotran ode
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

        field_states = field_states or []
        field_parameters = field_parameters or []
        monitored = monitored or []
        code_params = code_params or {}

        check_arg(ode, ODE)
        check_kwarg(field_states, "field_states", (list, tuple), itemtypes=str)
        check_kwarg(field_parameters, "field_parameters", (list, tuple), itemtypes=str)
        check_kwarg(monitored, "monitored", (list, tuple), itemtypes=str)
        check_kwarg(code_params, "code_params", dict)

        self.ode = ode
        self.field_states = field_states
        self.field_parameters = field_parameters
        self.monitored = monitored
        self._validate_input()

        params = GossCodeGeneratorParameters(**code_params)

        # Get a whole set of gotran code parameters and update with goss
        # specific options
        generation_params = parameters.generation.copy()
        generation_params.code.default_arguments = "st"
        generation_params.code.time.name = "time"

        generation_params.code.array.index_format = "[]"
        generation_params.code.array.index_offset = 0
        generation_params.code.array.flatten = True

        generation_params.code.parameters.representation = "named"
        generation_params.code.states.representation = params.state_repr.value
        generation_params.code.states.array_name = "states"

        generation_params.code.body.array_name = "body"
        generation_params.code.body.representation = params.body_repr.value
        generation_params.code.body.use_cse = params.use_cse
        generation_params.code.body.optimize_exprs = params.optimize_exprs.value
        generation_params.functions.jacobian.generate = params.generate_jacobian
        generation_params.functions.lu_factorization.generate = (
            params.generate_lu_factorization
        )
        generation_params.functions.forward_backward_subst.generate = (
            params.generate_forward_backward_subst
        )

        # Init base class
        super().__init__()

        # Store attributes
        self.params = generation_params
        self.file_form = _file_form.copy()
        self.class_form = _class_form.copy()

        self.name = ode.name.capitalize()

        if add_signature_to_name:
            signature = hashlib.sha1(
                (
                    ode.signature()
                    + repr(code_params)
                    + repr(field_states)
                    + repr(field_parameters)
                    + repr(monitored)
                    + gotran.__version__
                ).encode(),
            ).hexdigest()
            self.name += "_" + signature

        self._initialize_forms()

    def _initialize_forms(self):
        # Fill the forms with content
        self.file_form["MODELNAME"] = self.name.upper()
        self.class_form["ModelName"] = self.name
        self.class_form["num_states"] = self.ode.num_full_states
        self.class_form["num_parameters"] = self.ode.num_parameters
        self.class_form["num_field_states"] = len(self.field_states)
        self.class_form["num_field_parameters"] = len(self.field_parameters)
        self.class_form["num_monitored"] = len(self.monitored)
        self.class_form["monitored_evaluation_code"] = (
            _no_monitored_snippet.format(self.ode.name.capitalize()) + "\n"
        )
        self._code_generated = False

    def _validate_input(self):
        state_strs = [state.name for state in self.ode.full_states]
        for state_str in self.field_states:
            if state_str not in state_strs:
                error("{state_std} is not a state in the {self.ode} ODE")

        parameter_strs = [param.name for param in self.ode.parameters]
        for parameter_str in self.field_parameters:
            if parameter_str not in parameter_strs:
                error(f"{parameter_str} is not a parameter in the {self.ode} ODE")

        for expr_str in self.monitored:
            obj = self.ode.present_ode_objects.get(expr_str)
            if not isinstance(obj, Expression):
                error(f"{expr_str} is not an expression in the {self.ode} ODE")

    def class_code(self):
        """
        Generate the goss class code
        """

        if not self._code_generated:
            self._constructor_body()
            self._variable_init_and_declarations()
            self._eval_code()
            if len(self.monitored) > 0:
                self._monitored_code()

            self._jacobian_code()
            self._eval_linearized_code()
            self._eval_componentwise_code()

        self._code_generated = True

        return _class_template.format(**self.class_form)

    def file_code(self):
        """
        Generate the goss file code
        """
        self.file_form["CLASS_DECLARATION"] = self.class_code()
        return _file_template.format(**self.file_form)

    def _constructor_body(self):
        """
        Generate code snippets for constructor
        """

        ode = self.ode

        # State names
        state_names = [state.name for state in ode.full_states]
        body = ["", "// State names"]
        body.extend(
            '_state_names[{0}] = "{1}"'.format(i, name)
            for i, name in enumerate(state_names)
        )

        # Parameter names
        if self.class_form["num_parameters"] > 0:
            body.extend(["", "// Parameter names"])
            body.extend(
                '_parameter_names[{0}] = "{1}"'.format(i, param.name)
                for i, param in enumerate(ode.parameters)
            )

        # Field state names
        if self.class_form["num_field_states"] > 0:
            body.extend(["", "// Field state names"])
            body.extend(
                '_field_state_names[{0}] = "{1}"'.format(i, name)
                for i, name in enumerate(self.field_states)
            )

            body.extend(["", "// Field state indices"])
            for i, name in enumerate(self.field_states):
                body.append(
                    "_field_state_indices[{0}] = {1}".format(
                        i,
                        state_names.index(name),
                    ),
                )

        # Field parameter names
        if self.class_form["num_field_parameters"] > 0:
            body.extend(["", "// Field parameter names"])
            body.extend(
                '_field_parameter_names[{0}] = "{1}"'.format(i, name)
                for i, name in enumerate(self.field_parameters)
            )

        # Monitored names
        if self.class_form["num_monitored"] > 0:
            body.extend(["", "// Monitored names"])
            body.extend(
                '_monitored_names[{0}] = "{1}"'.format(i, monitored)
                for i, monitored in enumerate(self.monitored)
            )

        # Parameter to value map
        if self.class_form["num_parameters"] > 0:
            body.extend(["", "// Parameter to value map"])
            body.extend(
                '_param_to_value["{0}"] = &{1}'.format(param.name, param.name)
                for i, param in enumerate(ode.parameters)
            )

        # If ode is a DAE
        if ode.is_dae:
            body.extend(["", "// Set none differential states for DAE"])
            M = ode.mass_matrix
            for i in range(ode.num_full_states):
                if M[i, i].is_zero:
                    body.append("_differential_states[{0}] = 0".format(i))
            body.extend(["", "_is_dae = true"])

        body.append("")
        code = "\n".join(self.indent_and_split_lines(body, indent=3))

        self.class_form["constructor"] = code

    def _variable_init_and_declarations(self):
        """
        Generate code snippets for variable declarations, initialization and
        initial conditions
        """

        ode = self.ode

        parameter_declarations = []
        init = []
        init_copy = []

        # Parameter declaration and init
        if self.class_form["num_parameters"] > 0:
            parameter_declarations.extend(["", "// Parameters"])
            parameter_declarations.append(
                "double " + ", ".join(param.name for param in ode.parameters),
            )

            init.extend(
                "{0}({1})".format(param.name, param.init) for param in ode.parameters
            )
            init_copy.extend(
                "{0}(ode.{0})".format(param.name) for param in ode.parameters
            )

        # Parameter initialization
        init = [", ".join(init)]
        code = "\n".join(
            self.indent_and_split_lines(init, indent=3, no_line_ending=True),
        )
        if code:
            code = ",\n" + code
        self.class_form["variable_initialization"] = code

        init_copy = [", ".join(init_copy)]
        code = "\n".join(
            self.indent_and_split_lines(init_copy, indent=3, no_line_ending=True),
        )
        if code:
            code = ",\n" + code
        self.class_form["variable_initialization_copy_constructor"] = code

        # Parameter declaration
        code = "\n".join(self.indent_and_split_lines(parameter_declarations, indent=2))
        self.class_form["variable_declaration"] = code

        # Initial condition code
        ic_code = ["", "// Initial conditions"]
        ic_code.append("values->n = _num_states")
        ic_code.append("values->data.reset(new double[_num_states])")
        ic_code.extend(
            "values->data[{0}] = {1}".format(i, state.init)
            for i, state in enumerate(ode.full_states)
        )

        code = "\n".join(self.indent_and_split_lines(ic_code, indent=3))

        self.class_form["initial_condition_code"] = code

        # Field parameter setting
        set_field_parameters_code = []

        if self.class_form["num_field_parameters"] > 0:
            set_field_parameters_code.extend(["", "// Set field parameters"])
            set_field_parameters_code.extend(
                "{0} = field_params[{1}]".format(name, i)
                for i, name in enumerate(self.field_parameters)
            )

            set_field_parameters_code.append("")

            code = "\n".join(
                self.indent_and_split_lines(set_field_parameters_code, indent=3),
            )
            self.class_form["set_field_parameters_code"] = code

    def _eval_code(self):
        """
        Generate code for the eval method(s)
        """

        rhs = rhs_expressions(self.ode, result_name="values", params=self.params.code)
        self.class_form["eval_code"] = self.function_code(
            rhs,
            indent=3,
            include_signature=False,
        )

    def _monitored_code(self):
        """
        Generate code for the monitored method
        """
        monitored = monitored_expressions(
            self.ode,
            self.monitored,
            result_name="monitored",
            params=self.params.code,
        )
        self.class_form["monitored_evaluation_code"] = self.function_code(
            monitored,
            indent=3,
            include_signature=False,
        )

    def _eval_componentwise_code(self):
        body = ["", '//Timer timer_("Componentwise evaluation of rhs");', ""]
        body.extend(self.componentwise_code(self.ode, return_body_lines=True))
        body = self.wrap_body_with_function_prototype(
            body,
            "eval",
            "uint id, const double* states, double time",
            return_type="double",
            comment="Evaluate componenttwise rhs of the ODE",
        )
        code = "\n".join(self.indent_and_split_lines(body, indent=2))
        self.class_form["eval_componentwise_code"] = "\n" + code + "\n"

    def _eval_linearized_code(self):

        lin = linearized_derivatives(
            self.ode,
            function_name="linearized_eval",
            result_names=["linear", "rhs"],
            only_linear=False,
            include_rhs=True,
            nonlinear_last=True,
            params=self.params.code,
        )

        body = []
        for i in range(self.ode.num_full_states):
            if lin.linear_derivative_indices[i]:
                body.append("_linear_terms[{0}] = 1".format(i))

        code = "\n".join(self.indent_and_split_lines(body, indent=3))
        self.class_form["constructor_linear_terms"] = "\n" + code + "\n"

        body = ["", '//Timer timer_("Evaluation of linearized rhs");', ""]
        for line in self.function_code(lin, return_body_lines=True):
            if "Nonlinear linearized expressions" in line:
                body.append("// Return if only linear")
                body.append("if (only_linear)")
                body.append(["return"])
                body.append("")
            body.append(line)
        body = self.wrap_body_with_function_prototype(
            body,
            "linearized_eval",
            "const double* states, double time, "
            "double* linear, double* rhs, bool only_linear",
            comment="Evaluate the linearized rhs",
            const=True,
        )
        code = "\n".join(self.indent_and_split_lines(body, indent=2))
        self.class_form["linearized_eval_code"] = "\n" + code + "\n"

    def _jacobian_code(self):
        jac = None
        fact = None
        if self.params.functions.jacobian.generate:
            jac = jacobian_expressions(
                self.ode,
                result_name="jac",
                params=self.params.code,
            )

            body = ["", '//Timer timer_("Jacobian computation");', ""]
            body.extend(self.function_code(jac, return_body_lines=True))
            body = self.wrap_body_with_function_prototype(
                body,
                "compute_jacobian",
                "double* states, double time, double* jac",
                comment="Compute analytic jacobian",
            )
            code = "\n".join(self.indent_and_split_lines(body, indent=2))
            self.class_form["jacobian_code"] = "\n" + code + "\n"

        if self.params.functions.lu_factorization.generate:

            if jac is None:
                jac = jacobian_expressions(
                    self.ode,
                    result_name="jac",
                    params=self.params.code,
                )

            fact = factorized_jacobian_expressions(jac, params=self.params.code)
            body = ["", '//Timer timer_("Factorizing jacobian");', ""]
            body.extend(self.function_code(fact, return_body_lines=True))
            body.append("")
            body = self.wrap_body_with_function_prototype(
                body,
                "lu_factorize",
                "double* jac",
                comment="In place LU Factorize matrix (jacobian)",
                const=True,
            )
            code = "\n".join(self.indent_and_split_lines(body, indent=2))
            self.class_form["factorizing_code"] = "\n" + code + "\n"

        if self.params.functions.forward_backward_subst.generate:

            if jac is None:
                jac = jacobian_expressions(
                    self.ode,
                    result_name="jac",
                    params=self.params.code,
                )

            if fact is None:
                fact = factorized_jacobian_expressions(jac, params=self.params.code)

            fb_subst = forward_backward_subst_expressions(
                fact,
                residual_name="b",
                params=self.params.code,
            )
            body = ["", '//Timer timer_("Forward backward substitution");', ""]
            # body.extend(["// Copy b to dx", "for (unsigned int i=0; i<_num_states; i++)"])
            # body.append(["dx[i] = b[i]"])
            # body.append("")
            body.append("// In place forwards backward substitution")
            body.extend(self.function_code(fb_subst, return_body_lines=True))
            body.append("")
            body = self.wrap_body_with_function_prototype(
                body,
                "forward_backward_subst",
                "const double* jac, const double* b, double* dx",
                comment="Forward/Backward substitution of factoriesed matrix",
                const=True,
            )
            code = "\n".join(self.indent_and_split_lines(body, indent=2))
            self.class_form["fb_substitution_code"] = "\n" + code + "\n"
