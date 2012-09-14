#ifndef ExplicitEuler_h_IS_INCLUDED
#define ExplicitEuler_h_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>
#include <cstdlib>

namespace goss 
{

  // An Explicit Euler solver
  class ExplicitEuler : public ODESolver
  {
    
    public:
    
    // Default constructor
    ExplicitEuler() : ODESolver(0.0, 0.0), _dFdt(0) {};

    // Constructor
    ExplicitEuler(goss::ODE* ode_, double _ldt=-1.0);

    // Constructor
    ExplicitEuler(double _ldt);
    
    // Attach ODE to solver
    virtual void attach(goss::ODE* ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);

    // Destructor
    ~ExplicitEuler ();
    
  protected:

    // State derivative, allocated in attach(ode)
    double* _dFdt; 
    
  };
  
}
#endif
