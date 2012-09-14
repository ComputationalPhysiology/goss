#ifndef Sin_h_IS_INCLUDED
#define Sin_h_IS_INCLUDED

#include <goss/ODE.h>
#include <cmath>

namespace goss 
{

  class Sin : public ODE
  {
    public:
    double omega;
    Sin () : ODE(2, 0), omega(4*std::acos(0.0)) {}
    ~Sin() {}
    void eval(const double* y, double t, double* f_vals)
    {
      double y1 = y[0];
      double y2 = y[1];
      f_vals[0] =  omega*y2; 
      f_vals[1] = -omega*y1;
    }
    
    void getIC(goss::DoubleVector *res) const
    {
      res->n = system_size;
      res->data = new double[system_size];
      res->data[0] = 0.0; 
      res->data[1] = 1.0;//omega; 
    }
  };
}

#endif
