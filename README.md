[![CI-cpp](https://github.com/ComputationalPhysiology/goss/actions/workflows/cpp.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/cpp.yml)
[![CI-fenics](https://github.com/ComputationalPhysiology/goss/actions/workflows/fenics.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/fenics.yml)
[![github pages](https://github.com/ComputationalPhysiology/goss/actions/workflows/github-pages.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/github-pages.yml)
[![PyPI version](https://badge.fury.io/py/pygoss.svg)](https://badge.fury.io/py/pygoss)
[![coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/finsberg/a7290de789564f03eb6b1ee122fce423/raw/goss-badge.json)](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/finsberg/a7290de789564f03eb6b1ee122fce423/raw/goss-badge.json)

# `goss` - General ODE System Solver

`goss` is python wrapper around a C++ library for solving ordinary differential equations with a variety of different schemes.


Documentation is hosted at https://computationalphysiology.github.io/goss
Source code is found at https://github.com/ComputationalPhysiology/goss


## Motivation

The general idea is that you define your ODE in a [`gotran ode file`](https://github.com/ComputationalPhysiology/gotran) and hand the ode over to `goss`.

First define the ode in a gotran ODE file

```
# lorentz.ode
parameters(
sigma=10.0,
rho=28.0,
beta=8/3
)

# The values of the states represent the initial conditions
states(
x=0.0,
y=1.0,
z=1.05
)

dx_dt = sigma * (y - x)
dy_dt = x * (rho - z) - y
dz_dt = x * y - beta * z
```
You can now solve the ode as follows
```python
import numpy as np
import matplotlib.pyplot as plt
from gotran import load_ode
import goss

# Import the ode in gotran
lorentz = load_ode("lorentz.ode")
# Jit compile the code for the right hand side
ode = goss.ParameterizedODE(lorentz)
# Select a solver and instantiate the solver
solver = goss.solvers.RKF32(ode)
# Select the time steps you want to solve for
t = np.linspace(0, 100, 10001)
# Solve the system
u = solver.solve(t)

# Plot the solution
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

ax.plot(u[:, 0], u[:, 1], u[:, 2], lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Lorenz Attractor")
plt.show()
```
![_](https://raw.githubusercontent.com/ComputationalPhysiology/goss/main/docs/source/_static/lorentz.png)


For more examples, check out the [demo folder](https://github.com/ComputationalPhysiology/goss/tree/main/demo)

There is also a command line interface that can be used to list the available solvers, run the solver and generate the goss code

![_](https://raw.githubusercontent.com/ComputationalPhysiology/goss/main/docs/source/_static/cli.gif)



## Install

You can install goss with pip
```
python -m pip install pygoss
```
See [installation instructions](docs/install.md) for more options


## Known issues

- There is currently an issue on Apple Silicon with exceptions raised from by the jit compiled code which means the [one test](https://github.com/ComputationalPhysiology/goss/blob/main/tests/test_ode_bindings.py#L51) is not passing. An issue has been filed for this [here](https://github.com/wlav/cppyy/issues/68)

## Contributing

Contributions are very welcomed. To contribute please fork the repo, create a branch a submit a pull request. Before the pull request can be accepted it has to pass the test suit for the python and C++ code. Also note that we try to enforce an consistent coding style. To ensure that you follow the coding style you can install the pre-commit hook in the repo
```
python -m pip install pre-commit
pre-commit install
```
For every future commit, you will now run a set of tests that will make sure that you follow the coding style.

See the [contributing section](CONTRIBUTING.md) for more info.

## License
`goss` is licensed under the GNU LGPL, version 3 or (at your option) any later version.
