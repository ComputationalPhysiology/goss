from __future__ import annotations

from pathlib import Path
from typing import Any
from typing import Optional

import numpy as np
import rich_click as click
from gotran.model.loadmodel import load_ode
from modelparameters import utils

from . import solvers
from .codegeneration import GossCodeGenerator
from .codegeneration import GossCodeGeneratorParameters
from .ode import ParameterizedODE

try:
    import matplotlib.pyplot as plt

    has_mpl = True
except ImportError:
    plt = None
    has_mpl = False

try:
    from importlib.metadata import metadata
except ImportError:
    # python3.7 backport
    from importlib_metadata import metadata  # type: ignore

meta = metadata("pygoss")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]


class CodeParamsType(click.ParamType):
    name = "code-params"

    def convert(self, value, param, ctx):
        values = value.split(",")
        # Python use single quotes for the keys and we need double quotes
        # First convert ' -> $
        tmp = repr(dict(map(str.strip, v.split("=")) for v in values)).replace("'", "$")
        # Now convert ' > ""
        tmp = tmp.replace('"', "'")
        return GossCodeGeneratorParameters.parse_raw(tmp.replace("$", '"')).dict()


CodeParams = CodeParamsType()


@click.group()
@click.version_option(__version__, prog_name="goss")
def app():
    """
    goss - General ODE System Solver

    goss is a library for solving ODEs using the
    gotran ode format. It is written in C++ but
    contains python bindings.
    """
    pass


@click.command()
def list_solvers():
    """List available solvers and info about them"""
    from rich.console import Console
    from rich.table import Table

    table = Table(title="Goss solvers")

    table.add_column("Name", justify="right", style="cyan", no_wrap=True)
    table.add_column("Explicit/Implicit", style="magenta")
    table.add_column("Adaptive/Nonadaptive", style="green")

    for solver in solvers.GOSSSolvers._member_names_:
        exp_impl = (
            "Explicit"
            if solver in solvers.GOSSExplicitSolvers._member_names_
            else "Implicit"
        )
        adapt = (
            "Adaptive"
            if solver in solvers.GOSSAdaptiveSolvers._member_names_
            else "Nonadaptive"
        )
        table.add_row(solver, exp_impl, adapt)

    console = Console()
    console.print(table)


@click.command()
def code_params():
    """List available parameters to codegeneration"""
    GossCodeGeneratorParameters.print_defaults()


@click.command()
@click.argument(
    "filename",
    required=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
    ),
)
@click.option(
    "-o",
    "--output",
    default=None,
    type=click.Path(writable=True, resolve_path=True),
    help="Output file name",
)
@click.option(
    "--field-parameters",
    "-fp",
    type=str,
    help="Parameters which should be treated as field parameters",
    multiple=True,
)
@click.option(
    "--field-states",
    "-fs",
    type=str,
    help="States which should be treated as field states.",
    multiple=True,
)
@click.option(
    "--monitored",
    "-m",
    type=str,
    help="Monitored intermediates",
    multiple=True,
)
@click.option(
    "--list-timings",
    "-lt",
    type=bool,
    is_flag=True,
    default=False,
    help="Print timings of codegeneration.",
)
@click.option(
    "--code-params",
    "-cp",
    type=CodeParams,
    help="Parameters controling the code generation",
)
def gotran2goss(
    filename: Path,
    output: str = "",
    field_parameters: Optional[list[str]] = None,
    field_states: Optional[list[str]] = None,
    monitored: Optional[list[str]] = None,
    list_timings: bool = False,
    code_params: Optional[dict[str, Any]] = None,
):
    """
    Convert .ode file to goss file
    """
    field_parameters = field_parameters or []
    field_states = field_states or []
    monitored = monitored or []
    code_params = code_params or {}

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
        out = Path(filename).with_suffix(".h")

    with open(out, "w") as f:
        f.write(cgen.file_code())

    print(f"Output saved to {out}")

    if list_timings:
        utils.list_timings()


@click.command()
@click.argument(
    "filename",
    required=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
    ),
)
@click.option("--end-time", "-T", default=1000.0, type=float, help="End time")
@click.option(
    "--solver",
    type=click.Choice(solvers.GOSSSolvers._member_names_),
    default="ExplicitEuler",
    help="Which solver to use",
)
@click.option("--dt", "-dt", default=0.01, type=float, help="Time step")
@click.option(
    "--plot-y",
    "-y",
    type=str,
    help="States or monitored to plot on the y axis",
    multiple=True,
)
@click.option(
    "--plot-x",
    "-x",
    default="time",
    type=str,
    help="Values used for the x axis. Can be time and any valid plot_y variable.",
)
def run(
    filename: Path,
    end_time: float = 1000.0,
    solver: str = "ExplicitEuler",
    dt: float = 0.01,
    plot_y: Optional[list[str]] = None,
    plot_x: str = "time",
):
    """Solve an ODE"""
    if not has_mpl:
        click.echo("Matplotlib not installed - please install matplotlib")
        click.echo("python -m pip install matplotlib")
        click.Abort(ImportError)

    gotran_ode = load_ode(filename)

    if not plot_y:
        click.echo("Warning: ploy-y not specificed - assume you want to plot 'V'")
        plot_y = ["V"]  # Just assume you want to plot the membrane potential

    monitored_names = [
        expr.name for expr in gotran_ode.intermediates + gotran_ode.state_expressions
    ]
    monitored = [name for name in plot_y if name in monitored_names]
    if plot_x in monitored_names:
        monitored.append(plot_x)

    ode = ParameterizedODE(gotran_ode, monitored=monitored)

    cls = solvers.solver_mapper[solver]
    ode_solver = cls(ode)
    t = np.arange(0, end_time + dt, dt)
    y = ode_solver.solve(t)

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


app.add_command(list_solvers)
app.add_command(gotran2goss)
app.add_command(run)
app.add_command(code_params)
