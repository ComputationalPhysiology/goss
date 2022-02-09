#ifndef Arenstorf_h_IS_INCLUDED
#define Arenstorf_h_IS_INCLUDED

#include <memory>
#include <cmath>

#include <goss/ODE.h>

namespace goss {
  /*
    This is a system of ODEs describing a chemical reaction problem

    y1'=-0.04*y1+1.0e4*y2*y3              y1(0)=1
    y2'= 0.04*y1-1.0e4*y2*y3-3.0e7*y2*y2  y2(0)=0
    y3'= 3.0e7*y2*y2                      y3(0)=0

  */
  class Arenstorf : public ODE
  {
  private:
    double rval, rvalp;//parameters for the attractor

  public:

    Arenstorf (): ODE(4), rval(0.012277471), rvalp(1.0-rval)
    {}

    ~Arenstorf()
    {}

    std::shared_ptr<ODE> copy() const
    {
      return std::make_shared<Arenstorf>(*this);
    }

    virtual void eval(const double* y, double t, double* f_vals)
    {
      const double y0 = y[0];
      const double y1 = y[1];
      const double y2 = y[2];
      const double y3 = y[3];
      const double y1y1= y1*y1;
      double r1 = (y0 + rval)*(y0 + rval) + y1y1;
      double r2 = (y0 - rvalp)*(y0 - rvalp) + y1y1;
      r1 *= std::sqrt(r1);
      r2 *= std::sqrt(r2);
      f_vals[0] = y2;
      f_vals[1] = y3;
      f_vals[2] = y0 + 2*y3 - rvalp*(y0 + rval)/r1 - rval*(y0 - rvalp)/r2;
      f_vals[3] = y1 - 2*y2 - rvalp*y1/r1 - rval*y1/r2;
    }

    void get_ic(DoubleVector *res) const
    {
      res->n = _num_states;
      res->data.reset(new double[_num_states]);
      res->data[0] =  0.994;
      res->data[1] =  0.0;
      res->data[2] =  0.0;
      res->data[3] = -2.00158510637908252240537862224;
    }
  };
}

#endif
