#ifndef FHN_H_IS_INCLUDED
#define FHN_H_IS_INCLUDED
#include <memory>
#include <stdexcept>
#include <cmath>

#include "goss/ParameterizedODE.h"

namespace goss {

  // Implementation of gotran generated ODE
  class Fhn : public ParameterizedODE
  {
  public:

    // Constructor
    Fhn() : ParameterizedODE(2, 11, 0, 0, 0),
      a(0.13), b(0.013), c_1(0.26), c_2(0.1), c_3(1.0), stim_amplitude(50),
        stim_duration(1), stim_period(1000), stim_start(1), v_peak(40.0),
        v_rest(-85.0)

    {

      // State names
      _state_names[0] = "v";
      _state_names[1] = "s";

      // Parameter names
      _parameter_names[0] = "a";
      _parameter_names[1] = "b";
      _parameter_names[2] = "c_1";
      _parameter_names[3] = "c_2";
      _parameter_names[4] = "c_3";
      _parameter_names[5] = "stim_amplitude";
      _parameter_names[6] = "stim_duration";
      _parameter_names[7] = "stim_period";
      _parameter_names[8] = "stim_start";
      _parameter_names[9] = "v_peak";
      _parameter_names[10] = "v_rest";

      // Parameter to value map
      _param_to_value["a"] = &a;
      _param_to_value["b"] = &b;
      _param_to_value["c_1"] = &c_1;
      _param_to_value["c_2"] = &c_2;
      _param_to_value["c_3"] = &c_3;
      _param_to_value["stim_amplitude"] = &stim_amplitude;
      _param_to_value["stim_duration"] = &stim_duration;
      _param_to_value["stim_period"] = &stim_period;
      _param_to_value["stim_start"] = &stim_start;
      _param_to_value["v_peak"] = &v_peak;
      _param_to_value["v_rest"] = &v_rest;

      _linear_terms[1] = 1;

    }

    // Copy constructor
    Fhn(const Fhn& ode) : ParameterizedODE(ode),
      a(ode.a), b(ode.b), c_1(ode.c_1), c_2(ode.c_2), c_3(ode.c_3),
        stim_amplitude(ode.stim_amplitude), stim_duration(ode.stim_duration),
        stim_period(ode.stim_period), stim_start(ode.stim_start),
        v_peak(ode.v_peak), v_rest(ode.v_rest)
    {
      // Do nothing
    }

    // Evaluate rhs of the ODE
    void eval(const double* states, double time, double* values)
    {

      //Timer timer_("Evaluation of rhs");

      // Assign states
      const double v = states[0];
      const double s = states[1];

      // Expressions for the fhn component
      const double v_amp = v_peak - v_rest;
      const double v_th = a*v_amp + v_rest;
      const double I = -c_2*(v - v_rest)*s/v_amp + c_1*(v_peak - v)*(-v_th +
        v)*(v - v_rest)/(v_amp*v_amp);
      const double i_Stim = stim_amplitude*(1. - 1./(1. + std::exp(5.0*time -
        5.0*stim_start)))/(1. + std::exp(5.0*time - 5.0*stim_start -
        5.0*stim_duration));
      values[0] = I + i_Stim;
      values[1] = b*(v - c_3*s - v_rest);
    }


    // Evaluate the linearized rhs
    void linearized_eval(const double* states, double time, double* linear,
      double* rhs, bool only_linear) const
    {

      //Timer timer_("Evaluation of linearized rhs");


      // Assign states
      const double v = states[0];
      const double s = states[1];

      // Expressions for the fhn component
      const double v_amp = v_peak - v_rest;
      const double v_th = a*v_amp + v_rest;
      const double I = -c_2*(v - v_rest)*s/v_amp + c_1*(v_peak - v)*(-v_th +
        v)*(v - v_rest)/(v_amp*v_amp);
      const double i_Stim = stim_amplitude*(1. - 1./(1. + std::exp(5.0*time -
        5.0*stim_start)))/(1. + std::exp(5.0*time - 5.0*stim_start -
        5.0*stim_duration));
      rhs[0] = I + i_Stim;
      rhs[1] = b*(v - c_3*s - v_rest);
      linear[1] = -b*c_3;

      // Return if only linear
      if (only_linear)
      {
        return;
      }

      // Nonlinear linearized expressions
      const double dI_dv = -c_2*s/v_amp - c_1*(-v_th + v)*(v -
        v_rest)/(v_amp*v_amp) + c_1*(v_peak - v)*(-v_th + v)/(v_amp*v_amp) +
        c_1*(v_peak - v)*(v - v_rest)/(v_amp*v_amp);
      linear[0] = dI_dv;
    }

    // Evaluate componenttwise rhs of the ODE
    double eval(uint id, const double* states, double time)
    {

      //Timer timer_("Componentwise evaluation of rhs");

      // Return value
      double dy_comp[1] = {0.0};

      // What component?
      switch (id)
      {

        // Component 0 state v
        case 0:
        {

          // Assign states
          const double v = states[0];
          const double s = states[1];

          // Expressions for the fhn component
          const double v_amp = v_peak - v_rest;
          const double v_th = a*v_amp + v_rest;
          const double I = -c_2*(v - v_rest)*s/v_amp + c_1*(v_peak -
            v)*(-v_th + v)*(v - v_rest)/(v_amp*v_amp);
          const double i_Stim = stim_amplitude*(1. - 1./(1. +
            std::exp(5.0*time - 5.0*stim_start)))/(1. + std::exp(5.0*time -
            5.0*stim_start - 5.0*stim_duration));
          dy_comp[0] = I + i_Stim;
          break;
        }

        // Component 1 state s
        case 1:
        {

          // Assign states
          const double v = states[0];
          const double s = states[1];

          // Expressions for the fhn component
          dy_comp[0] = b*(v - c_3*s - v_rest);
          break;
        }

        // Default
        default:
        {
          throw std::runtime_error("Index out of bounds");
        }
      }

      // Return component
      return dy_comp[0];
    }

    // Get default initial conditions
    void get_ic(goss::DoubleVector *values) const
    {

      // Initial conditions
      values->n = _num_states;
      values->data.reset(new double[_num_states]);
      values->data[0] = -85.0;
      values->data[1] = 0.0;
    }

    // Return a copy of the ODE
    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<Fhn>(*this);
    }

    // Evaluate the monitored intermediates
    void eval_monitored(const double* states, double time, double* monitored) const
    {

      //Timer timer_("Evaluation of monitored.");

      // No monitored
      throw std::runtime_error("No monitored in the \'Fhn\' model.");

    }

    // Set all field parameters
    void set_field_parameters(const double* field_params)
    {

    }

  private:

    // Parameters
    double a, b, c_1, c_2, c_3, stim_amplitude, stim_duration, stim_period,
      stim_start, v_peak, v_rest;

  };

}

#endif
