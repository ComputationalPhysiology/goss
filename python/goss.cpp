#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_ODE(py::module &);

namespace goss {

PYBIND11_MODULE(goss, m) {
    // Optional docstring
    m.doc() = "goss library";
    
    init_ODE(m);
}
}