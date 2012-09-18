#ifndef RK4_h_IS_INCLUDED
#define RK4_h_IS_INCLUDED

#include "ODESolver.h"
#include "types.h"

namespace goss 
{
  
  // Explicit Runge Kutta solver of 4th order
  class RK4 : public ODESolver 
  {
  
  public:
  
    // Default constructor
    RK4();

    // Constructor
    RK4(double ldt);

    // Constructor
    RK4(goss::ODE* ode, double ldt=-1.0);

    // Destructor
    ~RK4();

    // Attach ODE to solver
    virtual void attach(goss::ODE* ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);

  protected: 
    
    // State derivative, allocated in attach(ode)
    double *k1, *k2, *k3, *k4, *tmp;

    // Perform a weighted addition of y and z
    inline void axpy(double* x, const double* y, double a, const double* z);
    
  };
}
//-----------------------------------------------------------------------------
inline void goss::RK4::axpy(double* x, const double* y, double a, const double* z)
{
  for (uint i = 0; i < ode_size(); ++i) 
    x[i] = y[i] + a*z[i];
}
//-----------------------------------------------------------------------------
#endif
