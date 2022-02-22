#ifndef OSCILATOR_H_IS_INCLUDED
#define OSCILATOR_H_IS_INCLUDED
#include <memory>
#include <stdexcept>
#include <cmath>

#include "goss/ParameterizedODE.h"

namespace goss {

  // Implementation of gotran generated ODE
  class Oscilator : public ParameterizedODE
  {
  public:

    // Constructor
    Oscilator() : ParameterizedODE(2, 2, 0, 0, 0),
      a(1.0), b(1.0)

    {

      // State names
      _state_names[0] = "x";
      _state_names[1] = "y";

      // Parameter names
      _parameter_names[0] = "a";
      _parameter_names[1] = "b";

      // Parameter to value map
      _param_to_value["a"] = &a;
      _param_to_value["b"] = &b;



    }

    // Copy constructor
    Oscilator(const Oscilator& ode) : ParameterizedODE(ode),
      a(ode.a), b(ode.b)
    {
      // Do nothing
    }

    // Evaluate rhs of the ODE
    void eval(const double* states, double time, double* values)
    {

      //Timer timer_("Evaluation of rhs");

      // Assign states
      const double x = states[0];
      const double y = states[1];

      // Expressions for the oscilator component
      values[0] = -a*y;
      values[1] = b*x;
    }


    // Evaluate the linearized rhs
    void linearized_eval(const double* states, double time, double* linear,
      double* rhs, bool only_linear) const
    {

      //Timer timer_("Evaluation of linearized rhs");


      // Assign states
      const double x = states[0];
      const double y = states[1];

      // Expressions for the oscilator component
      rhs[0] = -a*y;
      rhs[1] = b*x;

      // Return if only linear
      if (only_linear)
      {
        return;
      }

      // Nonlinear linearized expressions
      linear[0] = 0.;
      linear[1] = 0.;
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

        // Component 0 state x
        case 0:
        {

          // Assign states
          const double y = states[1];

          // Expressions for the oscilator component
          dy_comp[0] = -a*y;
          break;
        }

        // Component 1 state y
        case 1:
        {

          // Assign states
          const double x = states[0];

          // Expressions for the oscilator component
          dy_comp[0] = b*x;
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
      values->data[0] = 1.0;
      values->data[1] = 0.0;
    }

    // Return a copy of the ODE
    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<Oscilator>(*this);
    }

    // Evaluate the monitored intermediates
    void eval_monitored(const double* states, double time, double* monitored) const
    {

      //Timer timer_("Evaluation of monitored.");

      // No monitored
      throw std::runtime_error("No monitored in the \'Oscilator\' model.");

    }

    // Set all field parameters
    void set_field_parameters(const double* field_params)
    {

    }

  private:

    // Parameters
    double a, b;

  };

}

extern "C" DLL_EXPORT goss::ParameterizedODE * create_Oscilator()
{
  return new goss::Oscilator;
}


#endif
