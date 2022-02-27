from pathlib import Path
from unittest import mock

import goss
import gotran
import pytest
from click.testing import CliRunner

here = Path(__file__).absolute().parent
odefile = here.joinpath("fitzhughnagumo.ode")
runner = CliRunner()


@pytest.fixture
def output_file():
    output = here.joinpath("testfile.h")
    yield output
    output.unlink()


def test_gotran2goss(output_file):
    result = runner.invoke(
        goss.cli.app,
        ["gotran2goss", "--output", output_file.as_posix(), odefile.as_posix()],
    )
    assert result.exit_code == 0
    assert output_file.is_file()


def test_gotran2goss_fields_and_list_timings():
    result = runner.invoke(
        goss.cli.app,
        [
            "gotran2goss",
            "--list-timings",
            "-fp",
            "a",
            "-fp",
            "b",
            "-fs",
            "V",
            "-m",
            "I",
            odefile.as_posix(),
        ],
    )
    assert result.exit_code == 0
    output = odefile.with_suffix(".h")
    assert output.is_file()
    assert "num  : total time : mean time" in result.stdout
    output.unlink()


def test_gotran2goss_and_loadfile(output_file):
    result = runner.invoke(
        goss.cli.app,
        [
            "gotran2goss",
            "--output",
            output_file.as_posix(),
            odefile.as_posix(),
        ],
    )
    assert result.exit_code == 0
    gotran_ode = gotran.load_ode(odefile)
    ode = goss.ParameterizedODE(output_file, name=gotran_ode.name.capitalize())
    assert ode.num_states == gotran_ode.num_states


def test_gotran2goss_code_params():
    result = runner.invoke(
        goss.cli.app,
        [
            "gotran2goss",
            "-cp",
            "state_repr=named,generate_jacobian=True",
            odefile.as_posix(),
        ],
    )
    assert result.exit_code == 0
    output = odefile.with_suffix(".h")
    assert output.is_file()
    output.unlink()


def test_gossrun():
    with mock.patch("goss.cli.plt") as m:
        fig = mock.Mock()
        ax = mock.Mock()
        m.subplots.return_value = (fig, ax)
        result = runner.invoke(
            goss.cli.app,
            [
                "run",
                odefile.as_posix(),
                "--plot-y",
                "V",
                "--plot-y",
                "I",
                "--plot-x",
                "s",
            ],
        )
    assert result.exit_code == 0
    assert ax.plot.call_count == 2


def test_goss_list_solvers():
    result = runner.invoke(
        goss.cli.app,
        ["list-solvers"],
    )
    assert result.exit_code == 0
    assert "ExplicitEuler" in result.output


def test_goss_list_code_params():
    result = runner.invoke(
        goss.cli.app,
        ["code-params"],
    )
    assert result.exit_code == 0
    assert "GossCodeGeneratorParameters" in result.output
