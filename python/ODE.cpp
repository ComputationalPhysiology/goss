#include "../cpp/include/goss_bits/goss.h"

#include <pybind11/stl.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_ODE(py::module &m) {
    
    py::class_<goss::ODE>(m, "ODE")
    .def(py::init<uint>(), py::arg("num_states"))
    .def("num_states",
         py::overload_cast<>( &goss::ODE::num_states, py::const_));
}