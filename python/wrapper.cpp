#include <goss/goss.h>
#include <memory>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>


namespace py = pybind11;


void init_ODESolvers(py::module &m)
{
    py::class_<goss::ODESolver> solver_ODESolver(m, "ODESolver");

    solver_ODESolver.def("num_states", &goss::ODESolver::num_states)
            .def("is_adaptive", &goss::ODESolver::is_adaptive)
            .def_readwrite("ldt", &goss::ODESolver::ldt)
            .def("reset", &goss::ODESolver::reset)
            .def("get_internal_time_step", &goss::ODESolver::get_internal_time_step)
            .def("solve",
                 [](goss::ODESolver &self, py::array_t<double> y, py::array_t<double> y0,
                    py::array_t<double> t, const unsigned long num_timesteps, const int skip_n) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     py::buffer_info y0_info = y0.request();
                     auto y0_ptr = static_cast<double *>(y0_info.ptr);

                     py::buffer_info t_info = t.request();
                     auto t_ptr = static_cast<double *>(t_info.ptr);

                     self.solve(y_ptr, y0_ptr, t_ptr, num_timesteps, skip_n);
                 })
            .def("get_ode",
                 [](goss::ODESolver &self) {
                     return std::shared_ptr<const goss::ODE>(self.get_ode());
                 })
            .def("attach",
                 [](goss::ODESolver &self, std::shared_ptr<goss::ODE> ode) { self.attach(ode); });


    py::class_<goss::ExplicitEuler, goss::ODESolver> solver_ExplicitEuler(m, "ExplicitEuler");
    solver_ExplicitEuler.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("forward", [](goss::ExplicitEuler &self, const py::array_t<double> y, double t,
                               double interval) {
                py::buffer_info y_info = y.request();
                auto y_ptr = static_cast<double *>(y_info.ptr);

                self.forward(y_ptr, t, interval);
            });

    py::class_<goss::RL1, goss::ODESolver> solver_RL1(m, "RL1");
    solver_RL1.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("forward",
                 [](goss::RL1 &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 });

    py::class_<goss::GRL1, goss::ODESolver> solver_GRL1(m, "GRL1");
    solver_GRL1.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def_readwrite("delta", &goss::GRL1::delta)
            .def("forward",
                 [](goss::GRL1 &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 });


    // Implicit solvers

    py::class_<goss::ImplicitODESolver, goss::ODESolver> solver_ImplicitODESolver(
            m, "ImplicitODESolver");

    solver_ImplicitODESolver.def_readwrite("eta_0", &goss::ImplicitODESolver::eta_0)
            .def_readwrite("kappa", &goss::ImplicitODESolver::kappa)
            .def_readwrite("relative_tolerance", &goss::ImplicitODESolver::relative_tolerance)
            .def_readwrite("max_iterations", &goss::ImplicitODESolver::max_iterations)
            .def_readwrite("max_relative_previous_residual",
                           &goss::ImplicitODESolver::max_relative_previous_residual)
            .def_readwrite("always_recompute_jacobian",
                           &goss::ImplicitODESolver::always_recompute_jacobian)
            .def("num_jac_comp", &goss::ImplicitODESolver::num_jac_comp)
            .def("compute_factorized_jacobian",
                 [](goss::ImplicitODESolver &self, const py::array_t<double> y, double t, double dt,
                    double alpha) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.compute_factorized_jacobian(y_ptr, t, dt, alpha);
                 });

    py::class_<goss::ThetaSolver, goss::ImplicitODESolver> solver_ThetaSolver(m, "ThetaSolver");
    solver_ThetaSolver.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def_readwrite("num_refinements_without_always_recomputing_jacobian",
                           &goss::ThetaSolver::num_refinements_without_always_recomputing_jacobian)
            .def_readwrite("min_dt", &goss::ThetaSolver::min_dt)
            .def_readwrite("theta", &goss::ThetaSolver::theta)
            .def("forward", [](goss::ThetaSolver &self, const py::array_t<double> y, double t,
                               double interval) {
                py::buffer_info y_info = y.request();
                auto y_ptr = static_cast<double *>(y_info.ptr);

                self.forward(y_ptr, t, interval);
            });


    // Adaptive Implicit solvers

    py::class_<goss::AdaptiveImplicitSolver, goss::ImplicitODESolver> solver_AdaptiveImplicitSolver(
            m, "AdaptiveImplicitSolver");
    solver_AdaptiveImplicitSolver.def("get_atol", &goss::AdaptiveImplicitSolver::get_atol)
            .def("get_rtol", &goss::AdaptiveImplicitSolver::get_rtol)
            .def("get_iord", &goss::AdaptiveImplicitSolver::get_iord)
            .def("get_current_time", &goss::AdaptiveImplicitSolver::get_current_time)
            .def("get_current_time_step", &goss::AdaptiveImplicitSolver::get_current_time_step)
            .def("get_num_accepted", &goss::AdaptiveImplicitSolver::get_num_accepted)
            .def("get_num_rejected", &goss::AdaptiveImplicitSolver::get_num_rejected)
            .def("set_single_step_mode", [](goss::AdaptiveImplicitSolver &self,
                                            bool mode) { self.set_single_step_mode(mode); })
            .def("set_tol", [](goss::AdaptiveImplicitSolver &self, double atol,
                               double rtol = 1.0e-8) { self.set_tol(atol, rtol); })
            .def("set_iord",
                 [](goss::AdaptiveImplicitSolver &self, int iord) { self.set_iord(iord); });

    py::class_<goss::ESDIRK23a, goss::AdaptiveImplicitSolver> solver_ESDIRK23a(m, "ESDIRK23a");
    solver_ESDIRK23a.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def_readwrite("num_refinements_without_always_recomputing_jacobian",
                           &goss::ESDIRK23a::num_refinements_without_always_recomputing_jacobian)
            .def_readwrite("min_dt", &goss::ESDIRK23a::min_dt)
            .def_readonly("nfevals", &goss::ESDIRK23a::nfevals)
            .def_readonly("ndtsa", &goss::ESDIRK23a::ndtsa)
            .def_readonly("ndtsr", &goss::ESDIRK23a::ndtsr)
            .def("forward",
                 [](goss::ESDIRK23a &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 })
            .def("compute_ode_jacobian",
                 [](goss::ESDIRK23a &self, const py::array_t<double> y, double t) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.compute_ode_jacobian(y_ptr, t);
                 });
}


void init_ODE(py::module &m)
{

    m.def(
            "make_ode",
            [](std::uintptr_t e) {
                goss::ODE *p = reinterpret_cast<goss::ODE *>(e);
                return std::shared_ptr<const goss::ODE>(p);
            },
            "Create a goss::ODE object from a pointer integer, typically returned by a just-in-time compiler");

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

    py::class_<goss::ODE, PyODE, std::shared_ptr<goss::ODE>> ode(m, "ODE");


    // // py::class_<goss::ODE, PyODE>(m, "ODE")
    ode.def(py::init<goss::uint>())
            .def(py::init<const PyODE &>())
            .def("num_states", &goss::ODE::num_states)
            .def("is_dae", &goss::ODE::is_dae)
            .def("get_ic",
                 [](const goss::ODE &self) {
                     goss::DoubleVector values;
                     self.get_ic(&values);
                     return py::array_t<double>(values.n, values.data.get());
                 })
            .def("compute_jacobian",
                 [](goss::ODE &self, const py::array_t<double> states, double time,
                    const py::array_t<double> jac) {
                     py::buffer_info states_info = states.request();
                     auto states_ptr = static_cast<double *>(states_info.ptr);

                     py::buffer_info jac_info = jac.request();
                     auto jac_ptr = static_cast<double *>(jac_info.ptr);

                     self.compute_jacobian(states_ptr, time, jac_ptr);
                 })
            .def("linearized_eval",
                 [](goss::ODE &self, const py::array_t<double> states, double time,
                    const py::array_t<double> linearized, const py::array_t<double> rhs,
                    bool only_linear) {
                     py::buffer_info states_info = states.request();
                     auto states_ptr = static_cast<double *>(states_info.ptr);

                     py::buffer_info linearized_info = linearized.request();
                     auto linearized_ptr = static_cast<double *>(linearized_info.ptr);

                     py::buffer_info rhs_info = rhs.request();
                     auto rhs_ptr = static_cast<double *>(rhs_info.ptr);

                     self.linearized_eval(states_ptr, time, linearized_ptr, rhs_ptr, only_linear);
                 })
            .def("eval",
                 [](goss::ODE &self, const py::array_t<double> states, double time,
                    const py::array_t<double> values) {
                     py::buffer_info states_info = states.request();
                     auto states_ptr = static_cast<double *>(states_info.ptr);

                     py::buffer_info values_info = values.request();
                     auto values_ptr = static_cast<double *>(values_info.ptr);

                     self.eval(states_ptr, time, values_ptr);
                 })
            .def("eval",
                 [](goss::ODE &self, uint id, const py::array_t<double> states, double time) {
                     py::buffer_info states_info = states.request();
                     auto states_ptr = static_cast<double *>(states_info.ptr);

                     return self.eval(id, states_ptr, time);
                 });
}

PYBIND11_MODULE(_gosscpp, m)
{

    m.doc() = "This is a Python bindings of C++ goss Library";
    init_ODE(m);
    init_ODESolvers(m);
}
