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

try:
    import matplotlib.pyplot as plt

    has_mpl = True
except ImportError:
    has_mpl = False

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
    plot_y: List[str] = typer.Option(
        None,
        help="States or monitored to plot on the y axis.",
    ),
    plot_x: str = typer.Option(
        "time",
        help="Values used for the x axis. Can be time and any valid plot_y variable.",
    ),
):
    if not has_mpl:
        typer.echo("Matplotlib not installed - please install matplotlib")
        typer.echo("python -m pip install matplotlib")
        typer.Exit()

    gotran_ode = load_ode(filename)

    monitored_names = [
        expr.name for expr in gotran_ode.intermediates + gotran_ode.state_expressions
    ]
    monitored = [name for name in plot_y if name in monitored_names]
    if plot_x in monitored_names:
        monitored.append(plot_x)

    if len(plot_y) == 0:
        typer.echo("Warning: ploy-y not specificed - assume you want to plot 'V'")
        plot_y = ["V"]  # Just assume you want to plot the membrane potential

    ode = goss.ParameterizedODE(gotran_ode, monitored=monitored)

    cls = goss.solvers.solver_mapper[solver.name]
    ode_solver = cls(ode)
    y, t = ode_solver.solve(0, T, dt=dt)

    states = {}
    for name in plot_y:
        if name in gotran_ode.state_symbols:
            states[name] = gotran_ode.state_symbols.index(name)

    if len(monitored) > 0:
        m = ode.monitored_values(y, t)

    x = t
    if plot_x in gotran_ode.state_symbols:
        x = y[:, gotran_ode.state_symbols.index(plot_x)]
    if plot_x in monitored:
        x = m[:, monitored.index(plot_x)]

    for name in plot_y:
        fig, ax = plt.subplots()
        if name in gotran_ode.state_symbols:
            index = gotran_ode.state_symbols.index(name)
            ax.plot(x, y[:, index])
        elif name in monitored:
            index = monitored.index(name)
            ax.plot(x, m[:, index])

        ax.set_ylabel(name)
        ax.set_xlabel(plot_x)
    plt.show()
