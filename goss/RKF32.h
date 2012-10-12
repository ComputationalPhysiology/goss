#ifndef RKF32_h_IS_INCLUDED
#define RKF32_h_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "AdaptiveExplicitSolver.h"

namespace goss 
{

  // Adaptive and explicit RungeKutta Solver
  class RKF32 : public AdaptiveExplicitSolver
  {
  public:

    // Constructor
    RKF32();

    // Constructor
    RKF32(double _ldt);

    // Constructor
    RKF32 (boost::shared_ptr<ODE> ode, double ldt=-1.0);

    // Copy constructor
    RKF32(const RKF32& solver);

    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const { return boost::make_shared<RKF32>(*this); }

    // Constructor
    virtual ~RKF32();
    
    // Attach ODE
    virtual void attach(boost::shared_ptr<ODE> ode);

    // Reset solver
    virtual void reset();

    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);
    
    // FIXME: Where is this used!?
    // Store timestep and accepted timestep
    void log_data(double dt, bool accepted);
    
    // Return a vector of collected timesteps
    void dt_vector(DoubleVector *res);

    // Return a record of accepted 
    void accepted_vector(DoubleVector *res);

    // Counters for the number of right hand side evaluations (nfevals) and 
    // the number of accepted and rejected timesteps (ndtsa, ndtsr)
    long nfevals, ndtsa, ndtsr; 

  private: 

    // RK coefficients
    const double a21, a32;
    
    // RK weights
    const double b1, b2, b3, bh1, bh2, bh3, bh4;

    // Error weights
    const double d1, d2, d3, d4;

    // RK nodes
    const double c2, c3; 

    // System size in bytes
    ulong nbytes; 

    // State derivatives, allocated in attach(ode)
    boost::scoped_array<double> ki, k1, k2, k3, k4, yn, e;

    // Parameter for scalar or vector tolerance computing
    bool first;
  };
 
}
#endif
