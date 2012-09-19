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

    Sin () : ODE(2), omega(4*std::acos(0.0)) {}

    ODE* copy() const
    {
      return new Sin();
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
      res->n = _system_size;
      res->data.reset(new double[_system_size]);
      res->data[0] = 0.0; 
      res->data[1] = 1.0;//omega; 
    }
  };
}

#endif
