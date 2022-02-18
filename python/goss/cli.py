"""GOSS - General ODE solver
"""
from pathlib import Path
from typing import List
from typing import Optional

import goss
import typer
from goss.codegeneration import GossCodeGenerator
from gotran.model.loadmodel import load_ode
from modelparameters import utils

app = typer.Typer(no_args_is_help=True, help=__doc__)


def _list_solvers():
    from rich.console import Console
    from rich.table import Table

    table = Table(title="Goss solvers")

    table.add_column("Name", justify="right", style="cyan", no_wrap=True)
    table.add_column("Explicit/Implicit", style="magenta")
    table.add_column("Adaptive/Nonadaptive", style="green")

    for solver in goss.goss_solvers:
        exp_impl = "Explicit" if solver in goss.goss_explicit_solvers else "Implicit"
        adapt = "Adaptive" if solver in goss.goss_adaptive_solvers else "Nonadaptive"
        table.add_row(solver.__name__, exp_impl, adapt)

    console = Console()
    console.print(table)


@app.command("solvers", help="List available solvers and info about them")
def solvers(
    list_solvers: bool = typer.Option(
        False,
        "--list",
        "-l",
        help="List all available solvers",
    ),
):
    if list_solvers:
        _list_solvers()


@app.command("gotran2goss", help="Convert .ode file to goss file")
def gotran2goss(
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
