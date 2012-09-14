#ifndef VDP_h_IS_INCLUDED
#define VDP_h_IS_INCLUDED

#include <goss/ODE.h>

namespace goss 
{

  class VDP : public ODE
  {
  public:
    double mu;

    VDP() : ODE(2, 0), mu (10.0){}
    ~VDP() {}
    
    void eval(const double* y, double t, double* f_vals)
    {
      const double y1 = y[0];
      const double y2 = y[1];
      f_vals[0] = y2;
      f_vals[1] = mu*mu*((1-y1*y1)*y2-y1);
    }
    
    void getIC(goss::DoubleVector *res) const
    {
      res->n = system_size;
      res->data = new double[system_size];
      res->data[0] = 2.0;
      res->data[1] = 0.0;
    }

  };
}

#endif
