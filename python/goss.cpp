#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_goss(py::module &);

namespace mcl {

PYBIND11_MODULE(goss, m) {
    // Optional docstring
    m.doc() = "goss library";
    
    init_goss(m);
}
}