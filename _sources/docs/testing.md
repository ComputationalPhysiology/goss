# Testing

## Python

The tests for the python code can be found in the folder [tests](https://github.com/ComputationalPhysiology/goss/tree/main/tests) and run with `pytest`. To run the tests, please install the test dependencies
```
python -m pip install ".[test]"
```
and run the tests with
```
python -m pytest
```

## C++

The C++ source code for `goss` is found in the folder [cpp](https://github.com/ComputationalPhysiology/goss/tree/main/cpp). The C++ code also has a separate test suite that can be found in [cpp/tests](https://github.com/ComputationalPhysiology/goss/tree/main/cpp/tests). To run the tests you need to first build goss with the BUILD_TESTS flag enabled

```
cmake -B build-cpp -S cpp -DBUILD_TESTS=ON
cmake --build build-cpp
```
and now you can run the tests
```
cd build-cpp
ctest
```
