#ifndef van_der_Pol_h_IS_INCLUDED
#define van_der_Pol_h_IS_INCLUDED

#include <memory>
#include <goss/ODE.h>
#include <cmath>

namespace goss
{

  class VanDerPol : public ODE
  {
    public:

    VanDerPol() : ODE(2)
    {
      _differential_states[1] = 0;
      _is_dae = true;
    }

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<VanDerPol>(*this);
    }

    ~VanDerPol() {}

    void eval(const double* y, double t, double* f_vals)
    {
      const double y_ = y[0];
      const double z_ = y[1];
      f_vals[0] =  -z_;
      f_vals[1] = y_ - (z_*z_/3.-1.)*z_;
    }

    void compute_jacobian(double* states, double time, double* jac)
    {

      const double z = states[1];

      jac[0] = 0.0;
      jac[1] = -1.0;
      jac[2] = 1.0;
      jac[3] = -z*z + 1;
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      const double z = 2.0;
      res->data[0] = (z*z/3.-1.)*z;
      res->data[1] = z;
    }
  };
}

#endif
