from pathlib import Path

import goss
import pytest
from typer.testing import CliRunner

here = Path(__file__).absolute().parent
odefile = here.joinpath("tentusscher_2004_mcell.ode")
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


def test_gotran2goss_no_outout_and_list_timings():
    result = runner.invoke(
        goss.cli.app,
        ["gotran2goss", "--list-timings", odefile.as_posix()],
    )
    assert result.exit_code == 0
    output = odefile.with_suffix(".h")
    assert output.is_file()
    assert "num  : total time : mean time" in result.stdout
    output.unlink()
