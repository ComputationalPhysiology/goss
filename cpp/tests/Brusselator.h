#ifndef Brusselator_h_IS_INCLUDED
#define Brusselator_h_IS_INCLUDED

#include <memory>
#include <goss/ODE.h>

namespace goss
{
  class Brusselator : public ODE
  {
  public:

    Brusselator () : ODE(2) {}
    ~Brusselator() {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<Brusselator>(*this);
    }

    virtual void eval(const double* state, double t, double* f_vals)
    {
      const double y1 = state[0];
      const double y2 = state[1];
      f_vals[0] = 1+y1*y1*y2-4*y1;
      f_vals[1] = 3*y1-y1*y1*y2;
    }

    virtual void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] = 1.5;
      res->data[1] = 3.0;
    }
  };
}

#endif
