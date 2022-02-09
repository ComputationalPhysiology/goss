#ifndef PANFILOVCSE_H_IS_INCLUDED
#define PANFILOVCSE_H_IS_INCLUDED

#include <memory>
#include <stdexcept>
#include <cmath>

#include "goss/ParameterizedODE.h"

namespace goss {

  // Implementation of gotran generated ODE
  class PanfilovCSE : public ParameterizedODE
  {
  public:

    // Constructor
    PanfilovCSE() : ParameterizedODE(2, 3, 1, 0, 0),
      time_constant(1.0), v_peak(35.0), v_rest(-85.0)

    {

      // State names
      _state_names[0] = "e";
      _state_names[1] = "g";

      // Parameter names
      _parameter_names[0] = "time_constant";
      _parameter_names[1] = "v_peak";
      _parameter_names[2] = "v_rest";

      // Field state names
      _field_state_names[0] = "e";

      // Field state indices
      _field_state_indices[0] = 0;

      // Parameter to value map
      _param_to_value["time_constant"] = &time_constant;
      _param_to_value["v_peak"] = &v_peak;
      _param_to_value["v_rest"] = &v_rest;

    }

    // Evaluate rhs of the ODE
    void eval(const double* states, double t, double* values)
    {

      // Assign states
      const double e = states[0];
      const double g = states[1];

      // Common Sub Expressions
      const double cse_0 = -v_rest;
      const double cse_1 = e + cse_0;
      const double cse_2 = v_peak + cse_0;
      const double cse_3 = 1.0/cse_2;
      const double cse_4 = cse_1*cse_3;
      values[0] = -time_constant*cse_2*(g*cse_4 + 8.0*cse_4*(cse_4 -
        1.0)*(cse_4 - 0.1));
      values[1] = time_constant*(-8.0*e*(cse_4 - 1.1) - g)*(0.07*g/(e + 0.3)
        + 0.01);

    }

    // Get default initial conditions
    void get_ic(goss::DoubleVector *values) const
    {

      // Initial conditions
      values->n = _num_states;
      values->data.reset(new double[_num_states]);
      values->data[0] = 0.0;
      values->data[1] = 0.0;
    }

    // Return a copy of the ODE
    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<PanfilovCSE>(*this);
    }

    // Evaluate the monitored intermediates
    void eval_monitored(const double* states, double t, double* monitored) const
    {

      // No monitored
      throw std::runtime_error("No monitored in the \'Panfilovcse\' model.");

    }

    // Set all field parameters
    void set_field_parameters(const double* field_params)
    {

    }

  private:

    // Parameters
    double time_constant, v_peak, v_rest;

  };

}

#endif
