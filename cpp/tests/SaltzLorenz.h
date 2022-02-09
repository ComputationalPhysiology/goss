#ifndef SaltzLorenz_h_IS_INCLUDED
#define SaltzLorenz_h_IS_INCLUDED

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
  class SaltzLorenz : public ODE
  {
  protected:
    double rho, r, b;//parameters for the attractor

  public:

    SaltzLorenz() : ODE(3), rho(10.0), r(28.0), b(8.0/3.0)
    {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<SaltzLorenz>(*this);
    }

    ~SaltzLorenz() {}

    void eval(const double* y, double t, double* f_vals)
    {
      const double y1 = y[0];
      const double y2 = y[1];
      const double y3 = y[2];
      f_vals[0] = -rho*y1+rho*y2;
      f_vals[1] = -y1*y3+r*y1-y2;
      f_vals[2] =  y1*y2-b*y3;
    }

    void get_ic(goss::DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] = -8.0;
      res->data[1] =  8.0;
      res->data[2] = 27.0;
    }
  };
}

#endif
