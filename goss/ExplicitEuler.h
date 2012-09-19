#ifndef ExplicitEuler_h_IS_INCLUDED
#define ExplicitEuler_h_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include <math.h>
#include <cstdlib>

#include "ODESolver.h"

namespace goss 
{

  // An Explicit Euler solver
  class ExplicitEuler : public ODESolver
  {
    
    public:
    
    // Default constructor
    ExplicitEuler();

    // Constructor
    ExplicitEuler(double _ldt);
    
    // Constructor
    ExplicitEuler(ODE* ode_, double _ldt=-1.0);

    // Copy constructor
    ExplicitEuler(const ExplicitEuler& solver);

    // Attach ODE to solver
    virtual void attach(goss::ODE* ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);

    // Destructor
    ~ExplicitEuler ();
    
  protected:

    // State derivative, allocated in attach(ode)
    boost::scoped_array<double> _dFdt; 
    
  };
  
}
#endif
