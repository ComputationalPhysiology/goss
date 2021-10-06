#include <goss/goss.h>

#include <pybind11/stl.h>

#include <pybind11/pybind11.h>


namespace py = pybind11;

void init_ExplicitEuler(py::module &m) {
    
    py::class_<goss::ExplicitEuler>(m, "ExplicitEuler")
    .def(py::init<>());
}

PYBIND11_MODULE(cpp, m) {

    m.doc() = "This is a Python binding of C++ goss Library";

   
    init_ExplicitEuler(m);
}