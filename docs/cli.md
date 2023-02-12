# Command line interface

When you install `goss` you will also get access to a command line interface that can be used to solve ODEs and to list information. You can get the help text by typing
```
$ goss --help

 Usage: goss [OPTIONS] COMMAND [ARGS]...

 goss - General ODE System Solver
 goss is a library for solving ODEs using the gotran ode format. It is
 written in C++ but contains python bindings.

╭─ Options ──────────────────────────────────────────────────────────────╮
│ --version      Show the version and exit.                              │
│ --help         Show this message and exit.                             │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Commands ─────────────────────────────────────────────────────────────╮
│ code-params      List available parameters to codegeneration           │
│ gotran2goss      Convert .ode file to goss file                        │
│ list-solvers     List available solvers and info about them            │
│ run              Solve an ODE                                          │
╰────────────────────────────────────────────────────────────────────────╯
```

The CLI is built on [`click`](https://click.palletsprojects.com), [`rich`](https://github.com/Textualize/rich) and [`pydantic`](https://pydantic-docs.helpmanual.io)


![_](_static/cli.gif)


## Use command line interface for solving the Ten Tusscher model

To solve an ODE from the command line, we can use the `goss run` commend. You can see all the available options by typing `goss run --help`
```
$ goss run --help

 Usage: goss run [OPTIONS] FILENAME

 Solve an ODE

╭─ Options ──────────────────────────────────────────────────────────────╮
│ --end-time  -T   FLOAT                      End time                   │
│ --solver         [ExplicitEuler|RK2|RK4|RL  Which solver to use        │
│                  1|RL2|GRL1|GRL2|ImplicitE                             │
│                  uler|ThetaSolver|RKF32|ES                             │
│                  DIRK23a]                                              │
│ --dt        -dt  FLOAT                      Time step                  │
│ --plot-y    -y   TEXT                       States or monitored to     │
│                                             plot on the y axis         │
│ --plot-x    -x   TEXT                       Values used for the x      │
│                                             axis. Can be time and any  │
│                                             valid plot_y variable.     │
│ --help                                      Show this message and      │
│                                             exit.                      │
╰────────────────────────────────────────────────────────────────────────╯
```

Assume that we have the Ten Tusscher model defined in the file `tentusscher_panfilov_2006_M_cell.ode`. To get this file you could do the following


1. The model is available at [CellML](https://models.physiomeproject.org/workspace/tentusscher_panfilov_2006)

2. Download model and convert
   ```
   git clone https://models.physiomeproject.org/workspace/tentusscher_panfilov_2006
   cd tentusscher_panfilov_2006
   ```

3. Convert `.cellml` file to `.ode`
   ```
   gotran cellml2gotran ten_tusscher_model_2006_IK1Ko_M_units.cellml
   ```


Now, to solve the model we could for example do

```
goss run -T 1000 --solver GRL1 -dt 0.01 --plot-y V --plot-y Ca_i tentusscher_panfilov_2006_M_cell.ode
```
Here we have specified the we want to solve for 1000ms using the `GRL1` solver and an internal time step of 0.01. We have also specified that we want plot the voltage (`V`) and the intracellular calcium concentration (`Ca_i`)

![_](_static/tentusscher.png)
