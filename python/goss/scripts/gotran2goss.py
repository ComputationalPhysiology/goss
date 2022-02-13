from typing import List
from typing import Optional

__author__ = "Johan Hake (hake.dev@gmail.com)"
__date__ = "2012-09-20 -- 2015-01-15"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__ = "GNU LGPL Version 3.0 or later"

from pathlib import Path
import typer

from modelparameters import utils
from gotran.model.loadmodel import load_ode
from goss.codegeneration import GossCodeGenerator


app = typer.Typer(
    no_args_is_help=True,
    help="Convert .ode file to goss file",
)


@app.callback(invoke_without_command=True)
def main(
    filename: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Specify output file name",
    ),
    output: str = typer.Option("", help="Specify output file name"),
    field_parameters: Optional[List[str]] = typer.Option(
        None,
        help="Parameters which should be treated as field parameters.",
    ),
    field_states: Optional[List[str]] = typer.Option(
        None,
        help="States which should be treated as field states.",
    ),
    monitored: Optional[List[str]] = typer.Option(
        None,
        help="Monitored intermediates.",
    ),
    list_timings: bool = typer.Option(False, help="Print timings of codegeneration."),
):
    """
    Create a c header file from a gotran model
    """
    field_parameters = field_parameters or []
    field_states = field_states or []
    monitored = monitored or []

    # TODO: Make it possible to add this from the command line
    code_params = GossCodeGenerator.default_parameters()

    # Load Gotran model
    ode = load_ode(filename)

    # Code generators
    cgen = GossCodeGenerator(
        ode=ode,
        field_states=field_states,
        field_parameters=field_parameters,
        monitored=monitored,
        code_params=code_params,
    )

    if output:
        out = Path(output).with_suffix(".h")
    else:
        out = filename.with_suffix(".h")

    with open(out, "w") as f:
        f.write("#ifndef {}_H_IS_INCLUDED\n".format(ode.name.upper()))
        f.write("#define {}_H_IS_INCLUDED\n".format(ode.name.upper()))
        f.write("#include <memory>\n")
        f.write("#include <stdexcept>\n")
        f.write("#include <cmath>\n\n")

        f.write('#include "goss/ParameterizedODE.h"\n\n')
        f.write(cgen.class_code())
        f.write("\n#endif\n")

    if list_timings:
        utils.list_timings()
