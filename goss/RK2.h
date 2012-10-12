#ifndef RK2_H_IS_INCLUDED
#define RK2_H_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

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
    RK2(double ldt);

    // Constructor
    RK2(boost::shared_ptr<ODE> ode, double ldt=-1.0);

    // Copy constructor
    RK2(const RK2& solver);

    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const { return boost::make_shared<RK2>(*this); }

    // Destructor
    ~RK2();

    // Attach ODE to solver
    void attach(boost::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    virtual void forward(double* y, double t, double dt);

  protected: 

    // State derivative, allocated in attach(ode)
    boost::scoped_array<double> k1, tmp;
    
    // Perform a weighted addition of y and z
    inline void axpy(double* x, const double* y, double a, const double* z);

  };
}
//-----------------------------------------------------------------------------
inline void goss::RK2::axpy(double* x, const double* y, double a, const double* z)
{
  for (uint i = 0; i < num_states(); ++i) 
    x[i] = y[i] + a*z[i];
}
//-----------------------------------------------------------------------------
#endif
