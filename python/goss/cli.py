"""GOSS - General ODE solver
"""
import typer

from . import scripts

app = typer.Typer(no_args_is_help=True, help=__doc__)
app.add_typer(scripts.gotran2goss.app, name="gotran2goss")
