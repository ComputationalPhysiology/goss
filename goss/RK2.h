#ifndef RK2_H_IS_INCLUDED
#define RK2_H_IS_INCLUDED

#include "ODESolver.h"
#include "types.h"

namespace goss 
{

  // Explicit Runge Kutta solver of 2nd order
  class RK2 : public ODESolver 
  {
  
  public:
    
    // Default constructor
    RK2();

    // Constructor
    RK2(goss::ODE* ode_);
    
    // Destructor
    ~RK2();

    // Attach ODE to solver
    void attach(goss::ODE* ode_);

    // Step solver an interval in time forward
    virtual void forward(double* y, double t, double dt);

  protected: 

    // State derivative, allocated in attach(ode)
    double *k1, *tmp;
    
    // Perform a weighted addition of y and z
    inline void axpy(double* x, const double* y, double a, const double* z);

  };
}
//-----------------------------------------------------------------------------
inline void goss::RK2::axpy(double* x, const double* y, double a, const double* z)
{
  for (uint i = 0; i < ode_size(); ++i) 
    x[i] = y[i] + a*z[i];
}
//-----------------------------------------------------------------------------
#endif
