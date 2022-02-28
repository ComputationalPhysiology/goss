#include <goss/goss.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <vector>

namespace py = pybind11;
void init_ODESolvers(py::module &m)
{
    py::class_<goss::ODESolver, std::shared_ptr<goss::ODESolver>> solver_ODESolver(m, "ODESolver");

    solver_ODESolver.def("num_states", &goss::ODESolver::num_states)
            .def("is_adaptive", &goss::ODESolver::is_adaptive)
            .def("reset", &goss::ODESolver::reset)
            .def("get_internal_time_step", &goss::ODESolver::get_internal_time_step)
            .def("set_internal_time_step",
                 [](goss::ODESolver &self, double time_step) {
                     self.set_internal_time_step(time_step);
                 })
            .def("solve",
                 [](goss::ODESolver &self, py::array_t<double> y, py::array_t<double> y0,
                    py::array_t<double> t, const unsigned long num_timesteps) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     py::buffer_info y0_info = y0.request();
                     auto y0_ptr = static_cast<double *>(y0_info.ptr);

                     py::buffer_info t_info = t.request();
                     auto t_ptr = static_cast<double *>(t_info.ptr);

                     self.solve(y_ptr, y0_ptr, t_ptr, num_timesteps);
                 })
            .def("get_ode",
                 [](goss::ODESolver &self) {
                     return std::shared_ptr<const goss::ODE>(self.get_ode());
                 })
            .def("attach",
                 [](goss::ODESolver &self, std::shared_ptr<goss::ODE> ode) { self.attach(ode); });

    // Explicit Solvers

    py::class_<goss::ExplicitEuler, goss::ODESolver, std::shared_ptr<goss::ExplicitEuler>>
            solver_ExplicitEuler(m, "ExplicitEuler");
    solver_ExplicitEuler.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::ExplicitEuler::copy)
            .def("forward", [](goss::ExplicitEuler &self, const py::array_t<double> y, double t,
                               double interval) {
                py::buffer_info y_info = y.request();
                auto y_ptr = static_cast<double *>(y_info.ptr);

                self.forward(y_ptr, t, interval);
            });

    py::class_<goss::RL1, goss::ODESolver, std::shared_ptr<goss::RL1>> solver_RL1(m, "RL1");
    solver_RL1.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::RL1::copy)
            .def("forward",
                 [](goss::RL1 &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 });

    py::class_<goss::RL2, goss::RL1, std::shared_ptr<goss::RL2>> solver_RL2(m, "RL2");
    solver_RL2.def(py::init<>()).def(py::init<std::shared_ptr<goss::ODE>>());

    py::class_<goss::GRL1, goss::ODESolver, std::shared_ptr<goss::GRL1>> solver_GRL1(m, "GRL1");
    solver_GRL1.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def_readwrite("delta", &goss::GRL1::delta)
            .def("copy", &goss::GRL1::copy)
            .def("forward",
                 [](goss::GRL1 &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 });

    py::class_<goss::GRL2, goss::GRL1, std::shared_ptr<goss::GRL2>> solver_GRL2(m, "GRL2");
    solver_GRL2.def(py::init<>()).def(py::init<std::shared_ptr<goss::ODE>>());

    py::class_<goss::RK2, goss::ODESolver, std::shared_ptr<goss::RK2>> solver_RK2(m, "RK2");
    solver_RK2.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::RK2::copy)
            .def("forward",
                 [](goss::RK2 &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 });

    py::class_<goss::RK4, goss::ODESolver, std::shared_ptr<goss::RK4>> solver_RK4(m, "RK4");
    solver_RK4.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::RK4::copy)
            .def("forward",
                 [](goss::RK4 &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 });

    // Adaptive Explicit Solvers

    py::class_<goss::AdaptiveExplicitSolver, goss::ODESolver,
               std::shared_ptr<goss::AdaptiveExplicitSolver>>
            solver_AdaptiveExplicitSolver(m, "AdaptiveExplicitSolver");
    solver_AdaptiveExplicitSolver.def("get_atol", &goss::AdaptiveExplicitSolver::get_atol)
            .def("get_rtol", &goss::AdaptiveExplicitSolver::get_rtol)
            .def("get_iord", &goss::AdaptiveExplicitSolver::get_iord)
            .def("get_current_time", &goss::AdaptiveExplicitSolver::get_current_time)
            .def("get_current_time_step", &goss::AdaptiveExplicitSolver::get_current_time_step)
            .def("get_num_accepted", &goss::AdaptiveExplicitSolver::get_num_accepted)
            .def("get_num_rejected", &goss::AdaptiveExplicitSolver::get_num_rejected)
            .def("set_single_step_mode", [](goss::AdaptiveExplicitSolver &self,
                                            bool mode) { self.set_single_step_mode(mode); })
            .def("set_tol", [](goss::AdaptiveExplicitSolver &self, double atol,
                               double rtol = 1.0e-8) { self.set_tol(atol, rtol); })
            .def("set_iord",
                 [](goss::AdaptiveExplicitSolver &self, int iord) { self.set_iord(iord); });


    py::class_<goss::RKF32, goss::AdaptiveExplicitSolver, std::shared_ptr<goss::RKF32>>
            solver_RKF32(m, "RKF32");
    solver_RKF32.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::RKF32::copy)
            .def_readonly("nfevals", &goss::RKF32::nfevals)
            .def_readonly("ndtsa", &goss::RKF32::ndtsa)
            .def_readonly("ndtsr", &goss::RKF32::ndtsr)
            .def("forward",
                 [](goss::RKF32 &self, const py::array_t<double> y, double t, double interval) {
                     py::buffer_info y_info = y.request();
                     auto y_ptr = static_cast<double *>(y_info.ptr);

                     self.forward(y_ptr, t, interval);
                 });


    // Implicit solvers

    py::class_<goss::ImplicitODESolver, goss::ODESolver, std::shared_ptr<goss::ImplicitODESolver>>
            solver_ImplicitODESolver(m, "ImplicitODESolver");

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

    py::class_<goss::BasicImplicitEuler, goss::ImplicitODESolver,
               std::shared_ptr<goss::BasicImplicitEuler>>
            solver_BasicImplicitEuler(m, "BasicImplicitEuler");
    solver_BasicImplicitEuler.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::BasicImplicitEuler::copy)
            .def("forward", [](goss::BasicImplicitEuler &self, const py::array_t<double> y,
                               double t, double interval) {
                py::buffer_info y_info = y.request();
                auto y_ptr = static_cast<double *>(y_info.ptr);

                self.forward(y_ptr, t, interval);
            });

    py::class_<goss::ImplicitEuler, goss::ImplicitODESolver, std::shared_ptr<goss::ImplicitEuler>>
            solver_ImplicitEuler(m, "ImplicitEuler");
    solver_ImplicitEuler.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::ImplicitEuler::copy)
            .def_readwrite(
                    "num_refinements_without_always_recomputing_jacobian",
                    &goss::ImplicitEuler::num_refinements_without_always_recomputing_jacobian)
            .def_readwrite("min_dt", &goss::ImplicitEuler::min_dt)
            .def("forward", [](goss::ImplicitEuler &self, const py::array_t<double> y, double t,
                               double interval) {
                py::buffer_info y_info = y.request();
                auto y_ptr = static_cast<double *>(y_info.ptr);

                self.forward(y_ptr, t, interval);
            });

    py::class_<goss::ThetaSolver, goss::ImplicitODESolver, std::shared_ptr<goss::ThetaSolver>>
            solver_ThetaSolver(m, "ThetaSolver");
    solver_ThetaSolver.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::ThetaSolver::copy)
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

    py::class_<goss::AdaptiveImplicitSolver, goss::ImplicitODESolver,
               std::shared_ptr<goss::AdaptiveImplicitSolver>>
            solver_AdaptiveImplicitSolver(m, "AdaptiveImplicitSolver");
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

    py::class_<goss::ESDIRK23a, goss::AdaptiveImplicitSolver, std::shared_ptr<goss::ESDIRK23a>>
            solver_ESDIRK23a(m, "ESDIRK23a");
    solver_ESDIRK23a.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::ESDIRK23a::copy)
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

    py::class_<goss::ESDIRK4O32, goss::AdaptiveImplicitSolver, std::shared_ptr<goss::ESDIRK4O32>>
            solver_ESDIRK4O32(m, "ESDIRK4O32");
    solver_ESDIRK4O32.def(py::init<>())
            .def(py::init<std::shared_ptr<goss::ODE>>())
            .def("copy", &goss::ESDIRK4O32::copy)
            .def_readonly("nfevals", &goss::ESDIRK4O32::nfevals)
            .def_readonly("ndtsa", &goss::ESDIRK4O32::ndtsa)
            .def_readonly("ndtsr", &goss::ESDIRK4O32::ndtsr)
            .def("forward", [](goss::ESDIRK4O32 &self, const py::array_t<double> y, double t,
                               double interval) {
                py::buffer_info y_info = y.request();
                auto y_ptr = static_cast<double *>(y_info.ptr);

                self.forward(y_ptr, t, interval);
            });
}

void init_ODE(py::module &m)
{
    m.def(
            "make_ode",
            [](std::uintptr_t e) {
                goss::ParameterizedODE *p = reinterpret_cast<goss::ParameterizedODE *>(e);
                return std::shared_ptr<const goss::ParameterizedODE>(p);
            },
            "Create a goss::ODE object from a pointer integer, typically returned by "
            "a just-in-time compiler");

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
            .def("lu_factorize",
                 [](const goss::ODE &self, py::array_t<double> mat) {
                     py::buffer_info mat_info = mat.request();
                     auto mat_ptr = static_cast<double *>(mat_info.ptr);

                     self.lu_factorize(mat_ptr);
                 })
            .def("forward_backward_subst",
                 [](const goss::ODE &self, const py::array_t<double> mat,
                    const py::array_t<double> b, py::array_t<double> x) {
                     py::buffer_info mat_info = mat.request();
                     auto mat_ptr = static_cast<double *>(mat_info.ptr);

                     py::buffer_info b_info = b.request();
                     auto b_ptr = static_cast<double *>(b_info.ptr);

                     py::buffer_info x_info = x.request();
                     auto x_ptr = static_cast<double *>(x_info.ptr);

                     self.forward_backward_subst(mat_ptr, b_ptr, x_ptr);
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

    class PyParameterizedODE : public goss::ParameterizedODE
    {
      public:
        /* Inherit the constructors */
        using goss::ParameterizedODE::ParameterizedODE;

        void eval_monitored(const double *states, double t, double *monitored) const override
        {
            PYBIND11_OVERLOAD_PURE_NAME(void, goss::ParameterizedODE, "eval_monitored",
                                        eval_monitored, states, t, monitored);
        }

        void set_field_parameters(const double *field_params) override
        {
            PYBIND11_OVERLOAD_PURE_NAME(void, goss::ParameterizedODE, "set_field_parameters",
                                        set_field_parameters, field_params);
        }
        void eval(const double *states, double time, double *values) override
        {
            PYBIND11_OVERLOAD_PURE_NAME(void, goss::ODE, "eval", eval, states, time, values);
        }
        void get_ic(goss::DoubleVector *values) const override
        {
            PYBIND11_OVERLOAD_PURE_NAME(void, goss::ODE, "get_ic", get_ic, values);
        }
        std::shared_ptr<ODE> copy() const override
        {
            PYBIND11_OVERLOAD_PURE_NAME(std::shared_ptr<goss::ODE>, goss::ODE, "copy", copy, );
        }
    };

    py::class_<goss::ParameterizedODE, goss::ODE, PyParameterizedODE,
               std::shared_ptr<goss::ParameterizedODE>>
            parameterized_ode(m, "ParameterizedODE");
    parameterized_ode.def(py::init<const PyParameterizedODE &>())
            .def(py::init<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>())
            .def("num_field_states", &goss::ParameterizedODE::num_field_states)
            .def("num_parameters", &goss::ParameterizedODE::num_parameters)
            .def("num_field_parameters", &goss::ParameterizedODE::num_field_parameters)
            .def("num_monitored", &goss::ParameterizedODE::num_monitored)
            .def("eval_monitored",
                 [](goss::ParameterizedODE &self, const py::array_t<double> states, double t,
                    py::array_t<double> monitored) {
                     py::buffer_info states_info = states.request();
                     auto states_ptr = static_cast<double *>(states_info.ptr);

                     py::buffer_info monitored_info = monitored.request();
                     auto monitored_ptr = static_cast<double *>(monitored_info.ptr);

                     return self.eval_monitored(states_ptr, t, monitored_ptr);
                 })
            .def("monitored_values",
                 [](goss::ParameterizedODE &self, const py::array_t<double> states,
                    const py::array_t<double> t, py::array_t<double> monitored,
                    py::array_t<double> m) {
                     py::buffer_info states_info = states.request();
                     auto states_ptr = static_cast<double *>(states_info.ptr);

                     py::buffer_info t_info = t.request();
                     auto t_ptr = static_cast<double *>(t_info.ptr);

                     py::buffer_info monitored_info = monitored.request();
                     auto monitored_ptr = static_cast<double *>(monitored_info.ptr);

                     py::buffer_info m_info = m.request();
                     auto m_ptr = static_cast<double *>(m_info.ptr);

                     return self.monitored_values(states_ptr, t_ptr, monitored_ptr, m_ptr,
                                                  t_info.shape[0]);
                 })
            .def("set_field_parameters",
                 [](goss::ParameterizedODE &self, const py::array_t<double> field_params) {
                     py::buffer_info field_params_info = field_params.request();
                     auto field_params_ptr = static_cast<double *>(field_params_info.ptr);

                     self.set_field_parameters(field_params_ptr);
                 })
            .def("set_parameter", [](goss::ParameterizedODE &self, std::string name,
                                     double value) { self.set_parameter(name, value); })
            .def("get_parameter", [](const goss::ParameterizedODE &self,
                                     std::string name) { return self.get_parameter(name); })
            .def("get_state_names", &goss::ParameterizedODE::get_state_names)
            .def("get_field_state_names", &goss::ParameterizedODE::get_field_state_names)
            .def("get_parameter_names", &goss::ParameterizedODE::get_parameter_names)
            .def("get_field_parameter_names", &goss::ParameterizedODE::get_field_parameter_names)
            .def("get_field_state_indices", &goss::ParameterizedODE::get_field_state_indices)
            .def("get_monitored_names", &goss::ParameterizedODE::get_monitored_names);
}

void init_ODESystemSolver(py::module &m)
{
    py::class_<goss::ODESystemSolver> solver_ODESystemSolver(m, "ODESystemSolver");

    solver_ODESystemSolver
            .def(py::init<goss::uint, std::shared_ptr<goss::ODESolver>,
                          std::shared_ptr<goss::ParameterizedODE>>())
            .def("num_nodes", &goss::ODESystemSolver::num_nodes)
            .def("get_num_threads", &goss::ODESystemSolver::get_num_threads)
            .def("set_num_threads",
                 [](goss::ODESystemSolver &self, goss::uint num_threads) {
                     self.set_num_threads(num_threads);
                 })
            .def("reset_default", &goss::ODESystemSolver::reset_default)
            .def("ode", &goss::ODESystemSolver::ode)
            .def("solver", &goss::ODESystemSolver::solver)
            .def("solve",
                 [](goss::ODESystemSolver &self, py::array_t<double> field_states,
                    py::array_t<double> t, const unsigned long num_timesteps,
                    bool tangled_storage) {
                     py::buffer_info field_states_info = field_states.request();
                     auto field_states_ptr = static_cast<double *>(field_states_info.ptr);

                     py::buffer_info t_info = t.request();
                     auto t_ptr = static_cast<double *>(t_info.ptr);

                     self.solve(field_states_ptr, t_ptr, num_timesteps, tangled_storage);
                 })
            .def("forward", [](goss::ODESystemSolver &self, double t,
                               double interval) { self.forward(t, interval); })
            .def("get_field_states",
                 [](const goss::ODESystemSolver &self, py::array_t<double> system_field_states,
                    bool tangled_storage) {
                     py::buffer_info system_field_states_info = system_field_states.request();
                     auto system_field_states_ptr =
                             static_cast<double *>(system_field_states_info.ptr);

                     self.get_field_states(system_field_states_ptr, tangled_storage);
                 })
            .def("set_field_states",
                 [](goss::ODESystemSolver &self, const py::array_t<double> system_field_states,
                    bool tangled_storage) {
                     py::buffer_info system_field_states_info = system_field_states.request();
                     auto system_field_states_ptr =
                             static_cast<double *>(system_field_states_info.ptr);

                     self.set_field_states(system_field_states_ptr, tangled_storage);
                 })
            .def("get_field_parameters",
                 [](const goss::ODESystemSolver &self, py::array_t<double> system_field_parameters,
                    bool tangled_storage) {
                     py::buffer_info system_field_parameters_info =
                             system_field_parameters.request();
                     auto system_field_parameters_ptr =
                             static_cast<double *>(system_field_parameters_info.ptr);

                     self.get_field_parameters(system_field_parameters_ptr, tangled_storage);
                 })
            .def("set_field_parameters",
                 [](goss::ODESystemSolver &self, const py::array_t<double> system_field_parameters,
                    bool tangled_storage) {
                     py::buffer_info system_field_parameters_info =
                             system_field_parameters.request();
                     auto system_field_parameters_ptr =
                             static_cast<double *>(system_field_parameters_info.ptr);

                     self.set_field_parameters(system_field_parameters_ptr, tangled_storage);
                 })
            .def("get_field_state_components",
                 [](const goss::ODESystemSolver &self, py::array_t<double> component_field_states,
                    const goss::uint num_components, py::array_t<goss::uint> components,
                    bool tangled_storage) {
                     py::buffer_info component_field_states_info = component_field_states.request();
                     auto component_field_states_ptr =
                             static_cast<double *>(component_field_states_info.ptr);

                     py::buffer_info components_info = components.request();
                     auto components_ptr = static_cast<goss::uint *>(components_info.ptr);

                     self.get_field_state_components(component_field_states_ptr, num_components,
                                                     components_ptr, tangled_storage);
                 })
            .def("set_field_state_components",
                 [](goss::ODESystemSolver &self, const py::array_t<double> component_field_states,
                    const goss::uint num_components, py::array_t<goss::uint> components,
                    bool tangled_storage) {
                     py::buffer_info component_field_states_info = component_field_states.request();
                     auto component_field_states_ptr =
                             static_cast<double *>(component_field_states_info.ptr);

                     py::buffer_info components_info = components.request();
                     auto components_ptr = static_cast<goss::uint *>(components_info.ptr);

                     self.set_field_state_components(component_field_states_ptr, num_components,
                                                     components_ptr, tangled_storage);
                 })
            .def("states",
                 [](goss::ODESystemSolver &self) {
                     const double *states = self.states();
                     auto num_states = self.ode()->num_states();
                     auto num_nodes = self.num_nodes();
                     return py::array_t<double>(std::vector<ptrdiff_t>{num_nodes, num_states},
                                                states);
                 })
            .def("states_at_node",
                 [](goss::ODESystemSolver &self, goss::uint node) {
                     const double *states = self.states(node);
                     auto num_states = self.ode()->num_states();
                     return py::array_t<double>(num_states, states);
                 })
            .def("forward", [](goss::ODESystemSolver &self, double t, double interval) {
                self.forward(t, interval);
            });
}

PYBIND11_MODULE(_gosscpp, m)
{
    m.doc() = "This is a Python bindings of C++ goss Library";
    m.def("has_openmp", &has_openmp, "Method for checking if goss compiled with OpenMP support");
    init_ODE(m);
    init_ODESolvers(m);
    init_ODESystemSolver(m);
}
