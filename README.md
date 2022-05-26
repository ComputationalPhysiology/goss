[![CI-cpp](https://github.com/ComputationalPhysiology/goss/actions/workflows/cpp.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/cpp.yml)
[![CI-fenics](https://github.com/ComputationalPhysiology/goss/actions/workflows/fenics.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/fenics.yml)
[![github pages](https://github.com/ComputationalPhysiology/goss/actions/workflows/github-pages.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/github-pages.yml)
[![PyPI version](https://badge.fury.io/py/pygoss.svg)](https://badge.fury.io/py/pygoss)
[![codecov](https://codecov.io/gh/ComputationalPhysiology/goss/branch/main/graph/badge.svg?token=Z7DVGX7SUR)](https://codecov.io/gh/ComputationalPhysiology/goss)

# GOSS - General ODE System Solver

`goss` is a C++ library for solving ordinary differential equations.

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


## Documentation

Documentation is hosed at https://computationalphysiology.github.io/goss

## Install

### Installing dependencies
First, you need to install `boost`, i.e

(Debian)
```
apt-get install libboost-all-dev
```
(Fedora)
```
yum install boost-devel
```
Mac OSX
```
brew install boost
```

### Installing goss
Next you can install the python package.
```
python -m pip install pygoss
```

Alternatively you can clone the repo, cd into it at execute
```
python -m pip install .
```
or use
```
python -m pip install -e .
```
for an editable install.

## Testing

### Python

The tests for the python code can be found in the folder [tests](tests) and run with `pytest`. To run the tests, please install the test dependencies
```
python -m pip install ".[test]"
```
and run the tests with
```
python -m pytest
```

### C++

The C++ source code for `goss` is found in the folder [cpp](cpp). The C++ code also has a separate test suite that can be found in [cpp/tests](cpp/tests). To run the tests you need to first build goss with the BUILD_TESTS flag enabled

```
cmake -B build-cpp -S cpp -DBUILD_TESTS=ON
cmake --build build-cpp
```
and now you can run the tests
```
cd build-cpp
ctest
```

## Structure

The bindings between python and C++ uses [pybind11](https://pybind11.readthedocs.io/en/stable/) and all the bindings are found in the file [python/wrapper.cpp](python/wrapper.cpp).

The python package is built using [scikit-build](https://scikit-build.readthedocs.io/en/latest/index.html) which is a build system especially suited for python code with C++ extensions.

## Contributing

Contributions are very welcomed. To contribute please fork the repo, create a branch a submit a pull request. Before the pull request can be accepted it has to pass the test suit for the python and C++ code. Also note that we try to enforce an consistent coding style. To ensure that you follow the coding style you can install the pre-commit hook in the repo
```
python -m pip install pre-commit
pre-commit install
```
For every future commit, you will now run a set of tests that will make sure that you follow the coding style.
