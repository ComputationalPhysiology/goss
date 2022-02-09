#ifndef Robertson_h_IS_INCLUDED
#define Robertson_h_IS_INCLUDED

#include <memory>
#include <goss/ODE.h>

namespace goss
{

  /*
     This is a system of ODEs describing a chemical reaction problem

     y1'=-0.04*y1+1.0e4*y2*y3              y1(0)=1
     y2'= 0.04*y1-1.0e4*y2*y3-3.0e7*y2*y2  y2(0)=0
     y3'= 3.0e7*y2*y2                      y3(0)=0

   */
  class Robertson : public ODE
  {
  public:

    Robertson () : ODE(3) {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<Robertson>(*this);
    }

    ~Robertson() {}

    void eval(const double* y, double t, double* f_vals)
    {
      const double y1 = y[0];
      const double y2 = y[1];
      const double y3 = y[2];
      f_vals[0] = -0.04*y1+1.0e4*y2*y3;
      f_vals[1] =  0.04*y1-1.0e4*y2*y3-3.0e7*y2*y2;
      f_vals[2] =  3.0e7*y2*y2;
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] = 1.0;
      res->data[1] = 0.0;
      res->data[2] = 0.0;
    }
  };
}

#endif
