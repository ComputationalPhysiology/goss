#ifndef Sin_h_IS_INCLUDED
#define Sin_h_IS_INCLUDED

#include <memory>
#include <goss/ODE.h>
#include <cmath>

namespace goss
{

  class Sin : public ODE
  {
    public:
    double omega;

    Sin () : ODE(2), omega(4*std::acos(0.0)) {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<Sin>(*this);
    }

    ~Sin() {}

    void eval(const double* y, double t, double* f_vals)
    {
      double y1 = y[0];
      double y2 = y[1];
      f_vals[0] =  omega*y2;
      f_vals[1] = -omega*y1;
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] = 0.0;
      res->data[1] = 1.0;//omega;
    }
  };
}

#endif
