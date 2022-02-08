# GOSS - General ODE solver

`goss` is a C++ library for solver ordinary differential equations.

The general idea is that you define your ODE in a [`gotran ode file`](https://github.com/ComputationalPhysiology/gotran) and hand hand the ode over to `goss`.

## Install
To work with goss from python, you only need to install the python package. Clone the repo, cd into it at execute
```
python -m pip install .
```
or use 
```
python -m pip install -e .
```
for an editable install.

TODO: Add goss to pypi


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

## Contributing

The bindings between python and C++ uses [pybind11](https://pybind11.readthedocs.io/en/stable/) and all the bindings are found in the file [python/wrapper.cpp](python/wrapper.cpp).

The python package is built using [scikit-build](https://scikit-build.readthedocs.io/en/latest/index.html) which is a build system especially suited for python code with C++ extensions.
