[![CI-cpp](https://github.com/ComputationalPhysiology/goss/actions/workflows/cpp.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/cpp.yml)
[![CI-fenics](https://github.com/ComputationalPhysiology/goss/actions/workflows/fenics.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/fenics.yml)
[![CI-python](https://github.com/ComputationalPhysiology/goss/actions/workflows/python.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/python.yml)
[![github pages](https://github.com/ComputationalPhysiology/goss/actions/workflows/github-pages.yml/badge.svg)](https://github.com/ComputationalPhysiology/goss/actions/workflows/github-pages.yml)
[![PyPI version](https://badge.fury.io/py/pygoss.svg)](https://badge.fury.io/py/pygoss)
[![codecov](https://codecov.io/gh/ComputationalPhysiology/goss/branch/main/graph/badge.svg?token=Z7DVGX7SUR)](https://codecov.io/gh/ComputationalPhysiology/goss)

# GOSS - General ODE solver

`goss` is a C++ library for solver ordinary differential equations.

The general idea is that you define your ODE in a [`gotran ode file`](https://github.com/ComputationalPhysiology/gotran) and hand the ode over to `goss`.

## Documentation

Documentation is hosed at https://computationalphysiology.github.io/goss

## Install
To work with goss from python, you only need to install the python package.
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
