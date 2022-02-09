#ifndef BEELER_REUTER_1977_H_IS_INCLUDED
#define BEELER_REUTER_1977_H_IS_INCLUDED
#include <memory>
#include <stdexcept>
#include <cmath>

#include "goss/ParameterizedODE.h"

namespace goss {

  // Implementation of gotran generated ODE
  class Beeler_reuter_1977 : public ParameterizedODE
  {
  public:

    // Constructor
    Beeler_reuter_1977() : ParameterizedODE(8, 10, 0, 0, 0),
      E_Na(50), g_Na(0.04), g_Nac(3e-05), g_s(0.0005), IstimAmplitude(0.5),
        IstimEnd(50000), IstimPeriod(1000), IstimPulseDuration(1),
        IstimStart(10), C(0.01)

    {

      // State names
      _state_names[0] = "m";
      _state_names[1] = "h";
      _state_names[2] = "j";
      _state_names[3] = "Cai";
      _state_names[4] = "d";
      _state_names[5] = "f";
      _state_names[6] = "x1";
      _state_names[7] = "V";

      // Parameter names
      _parameter_names[0] = "E_Na";
      _parameter_names[1] = "g_Na";
      _parameter_names[2] = "g_Nac";
      _parameter_names[3] = "g_s";
      _parameter_names[4] = "IstimAmplitude";
      _parameter_names[5] = "IstimEnd";
      _parameter_names[6] = "IstimPeriod";
      _parameter_names[7] = "IstimPulseDuration";
      _parameter_names[8] = "IstimStart";
      _parameter_names[9] = "C";

      // Parameter to value map
      _param_to_value["E_Na"] = &E_Na;
      _param_to_value["g_Na"] = &g_Na;
      _param_to_value["g_Nac"] = &g_Nac;
      _param_to_value["g_s"] = &g_s;
      _param_to_value["IstimAmplitude"] = &IstimAmplitude;
      _param_to_value["IstimEnd"] = &IstimEnd;
      _param_to_value["IstimPeriod"] = &IstimPeriod;
      _param_to_value["IstimPulseDuration"] = &IstimPulseDuration;
      _param_to_value["IstimStart"] = &IstimStart;
      _param_to_value["C"] = &C;

      _linear_terms[0] = 1;
      _linear_terms[1] = 1;
      _linear_terms[2] = 1;
      _linear_terms[4] = 1;
      _linear_terms[5] = 1;
      _linear_terms[6] = 1;

    }

    // Copy constructor
    Beeler_reuter_1977(const Beeler_reuter_1977& ode) : ParameterizedODE(ode),
      E_Na(ode.E_Na), g_Na(ode.g_Na), g_Nac(ode.g_Nac), g_s(ode.g_s),
        IstimAmplitude(ode.IstimAmplitude), IstimEnd(ode.IstimEnd),
        IstimPeriod(ode.IstimPeriod),
        IstimPulseDuration(ode.IstimPulseDuration),
        IstimStart(ode.IstimStart), C(ode.C)
    {
      // Do nothing
    }

    // Evaluate rhs of the ODE
    void eval(const double* states, double time, double* values)
    {

      //Timer timer_("Evaluation of rhs");

      // Assign states
      const double m = states[0];
      const double h = states[1];
      const double j = states[2];
      const double Cai = states[3];
      const double d = states[4];
      const double f = states[5];
      const double x1 = states[6];
      const double V = states[7];

      // Expressions for the Sodium current component
      const double i_Na = (g_Nac + g_Na*(m*m*m)*h*j)*(-E_Na + V);

      // Expressions for the m gate component
      const double V_eff_0 = (std::fabs(47. + V) < 0.01 ? 0.01 : 47. + V);
      const double alpha_m = (-47. - V)/(-1. + std::exp(-0.1*V_eff_0));
      const double beta_m = 0.709552672748991*std::exp(-0.056*V);
      values[0] = -beta_m*m + (1. - m)*alpha_m;

      // Expressions for the h gate component
      const double alpha_h = 5.49796243870906e-10*std::exp(-0.25*V);
      const double beta_h = 1.7/(1. + 0.158025320889648*std::exp(-0.082*V));
      values[1] = (1. - h)*alpha_h - beta_h*h;

      // Expressions for the j gate component
      const double alpha_j = 1.86904730072229e-10*std::exp(-0.25*V)/(1. +
        1.67882752999566e-7*std::exp(-0.2*V));
      const double beta_j = 0.3/(1. + 0.0407622039783662*std::exp(-0.1*V));
      values[2] = (1. - j)*alpha_j - beta_j*j;

      // Expressions for the Slow inward current component
      const double E_s = -82.3 - 13.0287*std::log(0.001*Cai);
      const double i_s = g_s*(-E_s + V)*d*f;
      values[3] = 7.0e-6 - 0.07*Cai - 0.01*i_s;

      // Expressions for the d gate component
      const double alpha_d = 0.095*std::exp(1./20. - V/100.)/(1. +
        1.43328813856966*std::exp(-0.0719942404607631*V));
      const double beta_d = 0.07*std::exp(-44./59. - V/59.)/(1. +
        std::exp(11./5. + V/20.));
      values[4] = -beta_d*d + (1. - d)*alpha_d;

      // Expressions for the f gate component
      const double alpha_f = 0.012*std::exp(-28./125. - V/125.)/(1. +
        66.5465065250986*std::exp(0.149925037481259*V));
      const double beta_f = 0.0065*std::exp(-3./5. - V/50.)/(1. +
        std::exp(-6. - V/5.));
      values[5] = (1. - f)*alpha_f - beta_f*f;

      // Expressions for the Time dependent outward current component
      const double i_x1 = 0.00197277571153285*(-1. +
        21.7584023961971*std::exp(0.04*V))*std::exp(-0.04*V)*x1;

      // Expressions for the X1 gate component
      const double alpha_x1 =
        0.0311584109863426*std::exp(0.0826446280991736*V)/(1. +
        17.4117080633277*std::exp(0.0571428571428571*V));
      const double beta_x1 =
        0.000391646440562322*std::exp(-0.0599880023995201*V)/(1. +
        std::exp(-4./5. - V/25.));
      values[6] = (1. - x1)*alpha_x1 - beta_x1*x1;

      // Expressions for the Time independent outward current component
      const double V_eff_1 = (std::fabs(23. + V) < 0.01 ? 0.01 : 23. + V);
      const double i_K1 = 0.0035*(4.6 + 0.2*V)/(1. - std::exp(-0.04*V_eff_1))
        + 0.0035*(-4. +
        119.856400189588*std::exp(0.04*V))/(8.33113748768769*std::exp(0.04*V)
        + 69.4078518387552*std::exp(0.08*V));

      // Expressions for the Stimulus protocol component
      const double Istim = (time - IstimPeriod*std::floor((time -
        IstimStart)/IstimPeriod) - IstimStart <= IstimPulseDuration && time
        <= IstimEnd && time >= IstimStart ? IstimAmplitude : 0.);

      // Expressions for the Membrane component
      values[7] = (-i_K1 - i_s - i_Na - i_x1 + Istim)/C;
    }


    // Evaluate the linearized rhs
    void linearized_eval(const double* states, double time, double* linear,
      double* rhs, bool only_linear) const
    {

      //Timer timer_("Evaluation of linearized rhs");


      // Assign states
      const double m = states[0];
      const double h = states[1];
      const double j = states[2];
      const double Cai = states[3];
      const double d = states[4];
      const double f = states[5];
      const double x1 = states[6];
      const double V = states[7];

      // Expressions for the Sodium current component
      const double i_Na = (g_Nac + g_Na*(m*m*m)*h*j)*(-E_Na + V);

      // Expressions for the m gate component
      const double V_eff_0 = (std::fabs(47. + V) < 0.01 ? 0.01 : 47. + V);
      const double alpha_m = (-47. - V)/(-1. + std::exp(-0.1*V_eff_0));
      const double beta_m = 0.709552672748991*std::exp(-0.056*V);
      rhs[0] = -beta_m*m + (1. - m)*alpha_m;
      linear[0] = -beta_m - alpha_m;

      // Expressions for the h gate component
      const double alpha_h = 5.49796243870906e-10*std::exp(-0.25*V);
      const double beta_h = 1.7/(1. + 0.158025320889648*std::exp(-0.082*V));
      rhs[1] = (1. - h)*alpha_h - beta_h*h;
      linear[1] = -alpha_h - beta_h;

      // Expressions for the j gate component
      const double alpha_j = 1.86904730072229e-10*std::exp(-0.25*V)/(1. +
        1.67882752999566e-7*std::exp(-0.2*V));
      const double beta_j = 0.3/(1. + 0.0407622039783662*std::exp(-0.1*V));
      rhs[2] = (1. - j)*alpha_j - beta_j*j;
      linear[2] = -alpha_j - beta_j;

      // Expressions for the Slow inward current component
      const double E_s = -82.3 - 13.0287*std::log(0.001*Cai);
      const double i_s = g_s*(-E_s + V)*d*f;
      rhs[3] = 7.0e-6 - 0.07*Cai - 0.01*i_s;

      // Expressions for the d gate component
      const double alpha_d = 0.095*std::exp(1./20. - V/100.)/(1. +
        1.43328813856966*std::exp(-0.0719942404607631*V));
      const double beta_d = 0.07*std::exp(-44./59. - V/59.)/(1. +
        std::exp(11./5. + V/20.));
      rhs[4] = -beta_d*d + (1. - d)*alpha_d;
      linear[4] = -beta_d - alpha_d;

      // Expressions for the f gate component
      const double alpha_f = 0.012*std::exp(-28./125. - V/125.)/(1. +
        66.5465065250986*std::exp(0.149925037481259*V));
      const double beta_f = 0.0065*std::exp(-3./5. - V/50.)/(1. +
        std::exp(-6. - V/5.));
      rhs[5] = (1. - f)*alpha_f - beta_f*f;
      linear[5] = -beta_f - alpha_f;

      // Expressions for the Time dependent outward current component
      const double i_x1 = 0.00197277571153285*(-1. +
        21.7584023961971*std::exp(0.04*V))*std::exp(-0.04*V)*x1;

      // Expressions for the X1 gate component
      const double alpha_x1 =
        0.0311584109863426*std::exp(0.0826446280991736*V)/(1. +
        17.4117080633277*std::exp(0.0571428571428571*V));
      const double beta_x1 =
        0.000391646440562322*std::exp(-0.0599880023995201*V)/(1. +
        std::exp(-4./5. - V/25.));
      rhs[6] = (1. - x1)*alpha_x1 - beta_x1*x1;
      linear[6] = -beta_x1 - alpha_x1;

      // Expressions for the Time independent outward current component
      const double V_eff_1 = (std::fabs(23. + V) < 0.01 ? 0.01 : 23. + V);
      const double i_K1 = 0.0035*(4.6 + 0.2*V)/(1. - std::exp(-0.04*V_eff_1))
        + 0.0035*(-4. +
        119.856400189588*std::exp(0.04*V))/(8.33113748768769*std::exp(0.04*V)
        + 69.4078518387552*std::exp(0.08*V));

      // Expressions for the Stimulus protocol component
      const double Istim = (time - IstimPeriod*std::floor((time -
        IstimStart)/IstimPeriod) - IstimStart <= IstimPulseDuration && time
        <= IstimEnd && time >= IstimStart ? IstimAmplitude : 0.);

      // Expressions for the Membrane component
      rhs[7] = (-i_K1 - i_s - i_Na - i_x1 + Istim)/C;

      // Return if only linear
      if (only_linear)
      {
        return;
      }

      // Nonlinear linearized expressions
      const double di_s_dE_s = -g_s*d*f;
      const double dE_s_dCai = -13.0287/Cai;
      linear[3] = -0.07 - 0.01*dE_s_dCai*di_s_dE_s;
      const double dV_eff_1_dV = (std::fabs(23. + V) < 0.01 ? 0. : 1.);
      const double di_K1_dV = 0.0035*(-4. +
        119.856400189588*std::exp(0.04*V))*(-5.55262814710042*std::exp(0.08*V)
        -
        0.333245499507508*std::exp(0.04*V))/((8.33113748768769*std::exp(0.04*V)
        +
        69.4078518387552*std::exp(0.08*V))*(8.33113748768769*std::exp(0.04*V)
        + 69.4078518387552*std::exp(0.08*V))) + 0.0007/(1. -
        std::exp(-0.04*V_eff_1)) +
        0.0167798960265423*std::exp(0.04*V)/(8.33113748768769*std::exp(0.04*V)
        + 69.4078518387552*std::exp(0.08*V)) - 0.00014*(4.6 +
        0.2*V)*dV_eff_1_dV*std::exp(-0.04*V_eff_1)/((1. -
        std::exp(-0.04*V_eff_1))*(1. - std::exp(-0.04*V_eff_1)));
      const double di_Na_dV = g_Nac + g_Na*(m*m*m)*h*j;
      const double di_s_dV = g_s*d*f;
      const double di_x1_dV = 0.00171697791075903*x1 -
        7.89110284613141e-5*(-1. +
        21.7584023961971*std::exp(0.04*V))*std::exp(-0.04*V)*x1;
      const double di_K1_dV_eff_1 = -0.00014*(4.6 +
        0.2*V)*std::exp(-0.04*V_eff_1)/((1. - std::exp(-0.04*V_eff_1))*(1. -
        std::exp(-0.04*V_eff_1)));
      linear[7] = (-di_Na_dV - dV_eff_1_dV*di_K1_dV_eff_1 - di_K1_dV -
        di_x1_dV - di_s_dV)/C;
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

        // Component 0 state m
        case 0:
        {

          // Assign states
          const double m = states[0];
          const double V = states[7];

          // Expressions for the m gate component
          const double V_eff_0 = (std::fabs(47. + V) < 0.01 ? 0.01 : 47. + V);
          const double alpha_m = (-47. - V)/(-1. + std::exp(-0.1*V_eff_0));
          const double beta_m = 0.709552672748991*std::exp(-0.056*V);
          dy_comp[0] = -beta_m*m + (1. - m)*alpha_m;
          break;
        }

        // Component 1 state h
        case 1:
        {

          // Assign states
          const double h = states[1];
          const double V = states[7];

          // Expressions for the h gate component
          const double alpha_h = 5.49796243870906e-10*std::exp(-0.25*V);
          const double beta_h = 1.7/(1. +
            0.158025320889648*std::exp(-0.082*V));
          dy_comp[0] = (1. - h)*alpha_h - beta_h*h;
          break;
        }

        // Component 2 state j
        case 2:
        {

          // Assign states
          const double j = states[2];
          const double V = states[7];

          // Expressions for the j gate component
          const double alpha_j = 1.86904730072229e-10*std::exp(-0.25*V)/(1. +
            1.67882752999566e-7*std::exp(-0.2*V));
          const double beta_j = 0.3/(1. + 0.0407622039783662*std::exp(-0.1*V));
          dy_comp[0] = (1. - j)*alpha_j - beta_j*j;
          break;
        }

        // Component 3 state Cai
        case 3:
        {

          // Assign states
          const double Cai = states[3];
          const double d = states[4];
          const double f = states[5];
          const double V = states[7];

          // Expressions for the Slow inward current component
          const double E_s = -82.3 - 13.0287*std::log(0.001*Cai);
          const double i_s = g_s*(-E_s + V)*d*f;
          dy_comp[0] = 7.0e-6 - 0.07*Cai - 0.01*i_s;
          break;
        }

        // Component 4 state d
        case 4:
        {

          // Assign states
          const double d = states[4];
          const double V = states[7];

          // Expressions for the d gate component
          const double alpha_d = 0.095*std::exp(1./20. - V/100.)/(1. +
            1.43328813856966*std::exp(-0.0719942404607631*V));
          const double beta_d = 0.07*std::exp(-44./59. - V/59.)/(1. +
            std::exp(11./5. + V/20.));
          dy_comp[0] = -beta_d*d + (1. - d)*alpha_d;
          break;
        }

        // Component 5 state f
        case 5:
        {

          // Assign states
          const double f = states[5];
          const double V = states[7];

          // Expressions for the f gate component
          const double alpha_f = 0.012*std::exp(-28./125. - V/125.)/(1. +
            66.5465065250986*std::exp(0.149925037481259*V));
          const double beta_f = 0.0065*std::exp(-3./5. - V/50.)/(1. +
            std::exp(-6. - V/5.));
          dy_comp[0] = (1. - f)*alpha_f - beta_f*f;
          break;
        }

        // Component 6 state x1
        case 6:
        {

          // Assign states
          const double x1 = states[6];
          const double V = states[7];

          // Expressions for the X1 gate component
          const double alpha_x1 =
            0.0311584109863426*std::exp(0.0826446280991736*V)/(1. +
            17.4117080633277*std::exp(0.0571428571428571*V));
          const double beta_x1 =
            0.000391646440562322*std::exp(-0.0599880023995201*V)/(1. +
            std::exp(-4./5. - V/25.));
          dy_comp[0] = (1. - x1)*alpha_x1 - beta_x1*x1;
          break;
        }

        // Component 7 state V
        case 7:
        {

          // Assign states
          const double m = states[0];
          const double h = states[1];
          const double j = states[2];
          const double Cai = states[3];
          const double d = states[4];
          const double f = states[5];
          const double x1 = states[6];
          const double V = states[7];

          // Expressions for the Sodium current component
          const double i_Na = (g_Nac + g_Na*(m*m*m)*h*j)*(-E_Na + V);

          // Expressions for the Slow inward current component
          const double E_s = -82.3 - 13.0287*std::log(0.001*Cai);
          const double i_s = g_s*(-E_s + V)*d*f;

          // Expressions for the Time dependent outward current component
          const double i_x1 = 0.00197277571153285*(-1. +
            21.7584023961971*std::exp(0.04*V))*std::exp(-0.04*V)*x1;

          // Expressions for the Time independent outward current component
          const double V_eff_1 = (std::fabs(23. + V) < 0.01 ? 0.01 : 23. + V);
          const double i_K1 = 0.0035*(4.6 + 0.2*V)/(1. -
            std::exp(-0.04*V_eff_1)) + 0.0035*(-4. +
            119.856400189588*std::exp(0.04*V))/(8.33113748768769*std::exp(0.04*V)
            + 69.4078518387552*std::exp(0.08*V));

          // Expressions for the Stimulus protocol component
          const double Istim = (time - IstimPeriod*std::floor((time -
            IstimStart)/IstimPeriod) - IstimStart <= IstimPulseDuration &&
            time <= IstimEnd && time >= IstimStart ? IstimAmplitude : 0.);

          // Expressions for the Membrane component
          dy_comp[0] = (-i_K1 - i_s - i_Na - i_x1 + Istim)/C;
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
      values->data[0] = 0.011;
      values->data[1] = 0.988;
      values->data[2] = 0.975;
      values->data[3] = 0.0001;
      values->data[4] = 0.003;
      values->data[5] = 0.994;
      values->data[6] = 0.0001;
      values->data[7] = -84.624;
    }

    // Return a copy of the ODE
    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<Beeler_reuter_1977>(*this);
    }

    // Evaluate the monitored intermediates
    void eval_monitored(const double* states, double time, double* monitored) const
    {

      //Timer timer_("Evaluation of monitored.");

      // No monitored
      throw std::runtime_error("No monitored in the \'Beeler_reuter_1977\' model.");

    }

    // Set all field parameters
    void set_field_parameters(const double* field_params)
    {

    }

  private:

    // Parameters
    double E_Na, g_Na, g_Nac, g_s, IstimAmplitude, IstimEnd, IstimPeriod,
      IstimPulseDuration, IstimStart, C;

  };

}
#endif
