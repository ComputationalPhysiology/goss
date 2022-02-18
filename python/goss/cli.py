"""GOSS - General ODE solver
"""
from pathlib import Path
from typing import List
from typing import Optional

import goss.solvers
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

    for solver in goss.solvers.GOSSSolvers._member_names_:
        exp_impl = (
            "Explicit"
            if solver in goss.solvers.GOSSExplicitSolvers._member_names_
            else "Implicit"
        )
        adapt = (
            "Adaptive"
            if solver in goss.solvers.GOSSIAdaptiveSolvers._member_names_
            else "Nonadaptive"
        )
        table.add_row(solver, exp_impl, adapt)

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
        for line in [
            f"#ifndef {ode.name.upper()}_H_IS_INCLUDED\n",
            f"#define {ode.name.upper()}_H_IS_INCLUDED\n",
            "#include <memory>\n",
            "#include <stdexcept>\n",
            "#include <cmath>\n\n",
            '#include "goss/ParameterizedODE.h"\n\n',
            cgen.class_code(),
            "\n#endif\n",
        ]:
            f.write(line)

    if list_timings:
        utils.list_timings()


@app.command("run", help="Solve an ODE")
def run(
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
    T: float = typer.Option(1000.0, "-T", help="End time"),
    solver: goss.solvers.GOSSSolvers = typer.Option(
        "ExplicitEuler",
        help="Which solver to use",
    ),
    dt: float = typer.Option(0.01, "-dt", help="Time step"),
):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        typer.echo("Please install matplotlib - pip install matplotlib")
        typer.Exit()

    gotran_ode = load_ode(filename)
    ode = goss.ODE(gotran_ode)
    cls = goss.solvers.solver_mapper[solver.name]
    ode_solver = cls(ode)
    y, t = ode_solver.solve(0, T, dt=dt)

    V_index = gotran_ode.state_symbols.index("V")

    fig, ax = plt.subplots()
    ax.plot(t, y[:, V_index])
    ax.set_title("V")
    plt.show()
