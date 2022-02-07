#include <boost/shared_ptr.hpp>
#include <goss/goss.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

void init_ExplicitEuler(py::module &m)
{

    py::class_<goss::ExplicitEuler>(m, "ExplicitEuler").def(py::init<>());
}


void init_ODE(py::module &m)
{

    class PyODE : public goss::ODE
    {
      public:
        /* Inherit the constructors */
        using goss::ODE::ODE;

        /* Trampoline (need one for each virtual function) */
        void eval(const double *states, double time, double *values) override
        {
            PYBIND11_OVERLOAD_PURE_NAME(void, goss::ODE, "eval", eval, states, time, values);
        }
        double eval(uint idx, const double *states, double time) override
        {
            PYBIND11_OVERLOAD_NAME(double, goss::ODE, "eval", eval, idx, states, time);
        }
        void linearized_eval(const double *states, double time, double *linearized, double *rhs,
                             bool only_linear) const override
        {
            PYBIND11_OVERLOAD_NAME(void, goss::ODE, "linearized_eval", linearized_eval, states,
                                   time, linearized, rhs, only_linear);
        }
        void get_ic(goss::DoubleVector *values) const override
        {
            PYBIND11_OVERLOAD_PURE_NAME(void, goss::ODE, "get_ic", get_ic, values);
        }
        std::shared_ptr<ODE> copy() const override
        {
            PYBIND11_OVERLOAD_PURE_NAME(std::shared_ptr<goss::ODE>, goss::ODE, "copy", copy, );
        }
        void compute_jacobian(double *states, double time, double *jac) override
        {
            PYBIND11_OVERLOAD_NAME(void, goss::ODE, "compute_jacobian", compute_jacobian, states,
                                   time, jac);
        }
        void lu_factorize(double *mat) const override
        {
            PYBIND11_OVERLOAD_NAME(void, goss::ODE, "lu_factorize", lu_factorize, mat);
        }
        void forward_backward_subst(const double *mat, const double *b, double *x) const override
        {
            PYBIND11_OVERLOAD_NAME(void, goss::ODE, "forward_backward_subst",
                                   forward_backward_subst, mat, b, x);
        }
    };

    py::class_<goss::ODE, PyODE, std::shared_ptr<goss::ODE>>ode(m, "ODE");


    // // py::class_<goss::ODE, PyODE>(m, "ODE")
    ode
      .def(py::init<goss::uint>())
      .def(py::init<const PyODE&>())
      .def("num_states", &goss::ODE::num_states);

}

PYBIND11_MODULE(_gosscpp, m)
{

    m.doc() = "This is a Python bindings of C++ goss Library";


    init_ExplicitEuler(m);
    init_ODE(m);
}