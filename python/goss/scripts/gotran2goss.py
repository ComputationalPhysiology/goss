#!/usr/bin/env python
__author__ = "Johan Hake (hake.dev@gmail.com)"
__date__ = "2012-09-20 -- 2015-01-15"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__ = "GNU LGPL Version 3.0 or later"

from pathlib import Path
import typer

from gotran.model.loadmodel import load_ode

# from gotran import list_timings
from goss.codegeneration import GossCodeGenerator


app = typer.Typer(
    no_args_is_help=True,
    help="Convert .ode file to goss file",
)
_code_params = GossCodeGenerator.default_parameters()


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
    )
):
    """
    Create a c header file from a gotran model
    """
    params = {
        "field_states": [],
        "field_parameters": [],
        "monitored": [],
        "output": "",
        "list_timings": False,
        "code": GossCodeGenerator.default_parameters(),
    }

    # Load Gotran model
    ode = load_ode(filename)

    # Code generators
    cgen = GossCodeGenerator(
        ode,
        params["field_states"],
        params["field_parameters"],
        params["monitored"],
        params["code"],
    )
    output = params["output"]

    if output:
        output = Path(output).with_suffix(".h")
    else:
        output = filename.with_suffix(".h")

    f = open(output, "w")
    f.write("#ifndef {}_H_IS_INCLUDED\n".format(ode.name.upper()))
    f.write("#define {}_H_IS_INCLUDED\n".format(ode.name.upper()))
    f.write("#include <boost/make_shared.hpp>\n")
    f.write("#include <stdexcept>\n")
    f.write("#include <cmath>\n\n")

    f.write('#include "goss/ParameterizedODE.h"\n\n')
    f.write(cgen.class_code())
    f.write("\n#endif\n")

    if params["list_timings"]:
        list_timings()


if __name__ == "__main__":

    params = ParameterDict(
        output=Param("", description="Specify output file name"),
        field_parameters=Param(
            [], description="Parameters which should be " "treated as field parameters."
        ),
        field_states=Param(
            [], description="States which should be " "treated as field states."
        ),
        monitored=Param([], description="Monitored intermediates."),
        list_timings=Param(False, description="Print timings of codegeneration."),
        code=code_params,
    )

    params.parse_args(usage="usage: %prog FILE [options]")  # sys.argv[2:])

    if len(sys.argv) < 2:
        error("Expected a single gotran file argument")

    if not os.path.isfile(sys.argv[1]):
        error("Expected the argument to be a file")

    file_name = sys.argv[1]
    main(file_name, params)
