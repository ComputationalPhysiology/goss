#ifndef ESDIRK4O32_h_IS_INCLUDED
#define ESDIRK4O32_h_IS_INCLUDED

#include "AdaptiveImplicitSolver.h"

namespace goss 
{

  // Explicit Singly Diagonally Implicit Runge-Kutta solver
  class ESDIRK4O32: public AdaptiveImplicitSolver
  {
  
    public:
    
    // Default constructor
    ESDIRK4O32();
    
    // Constructor
    ESDIRK4O32 (ODE* ode, double ldt=-1.0);

    // Constructor
    ESDIRK4O32(double _ldt);

    // Initilize data
    void init()
    {
      _iord = 3;
    }

    // Attach ODE
    virtual void attach(ODE* ode_);
    
    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);

    // Destructor
    ~ESDIRK4O32 ();

    // Counters for the number of right hand side evaluations (nfevals) and 
    // the number of accepted and rejected timesteps (ndtsa, ndtsr)
    long nfevals, ndtsa, ndtsr; 

  private:
    
    // Help variable
    double gamma;
    
    // RK coefficients
    double a21, a22, a31, a32, a33, a41, a42, a43, a44;

    // RK weights
    double b1, b2, b3, b4, bh1, bh2, bh3;
    
    // RK coefficients
    double c2, c3, c4;

    // State derivatives, allocated in attach(ode)
    double *z1, *z2, *z3, *z4, *yn, *yh, *swap, *ret_ptr; 

  };

}
#endif
