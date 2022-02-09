#ifndef NonLinOscillator_h_IS_INCLUDED
#define NonLinOscillator_h_IS_INCLUDED

#include <memory>
#include <cmath>

#include <goss/ODE.h>

namespace goss
{

  class NonLinOscillator : public ODE
  {
  public:
    double omega;

    NonLinOscillator() : ODE(2), omega(4*std::acos(0.0)) {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<NonLinOscillator>(*this);
    }

    ~NonLinOscillator() {}

    void eval(const double* y, double t, double* f_vals)
    {
      const double y1 = y[0];
      const double y2 = y[1];
      f_vals[0] = sin(y2)*(y1-6.0)/6.0;
      f_vals[1] = y1+cos(y2)/6.0;
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] = 0.0;
      res->data[1] = acos(-0.8);//omega;
    }
  };
}

#endif
