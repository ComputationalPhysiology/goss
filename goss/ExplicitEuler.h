#ifndef ExplicitEuler_h_IS_INCLUDED
#define ExplicitEuler_h_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
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
    ExplicitEuler(boost::shared_ptr<ODE> ode, double _ldt=-1.0);

    // Copy constructor
    ExplicitEuler(const ExplicitEuler& solver);

    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const 
    { return boost::make_shared<ExplicitEuler>(*this); }

    // Attach ODE to solver
    virtual void attach(boost::shared_ptr<ODE> ode);

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
