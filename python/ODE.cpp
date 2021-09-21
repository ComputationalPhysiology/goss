#include "../cpp/include/goss_bits/goss.h"

#include <pybind11/stl.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_ExplicitEuler(py::module &m) {
    
    py::class_<goss::ExplicitEuler>(m, "ExplicitEuler")
    .def(py::init<>());
}