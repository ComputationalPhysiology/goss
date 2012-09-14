#ifndef ODESolver_h_IS_INCLUDED
#define ODESolver_h_IS_INCLUDED

#include <iostream>

#include "ODE.h"
#include "types.h"

namespace goss 
{

  class ODESolver {

  public:

    // FIXME: Remove dt as variable...
    // Constructor
    ODESolver (goss::ODE* ode, double ldt=0.0, double dt=0.0)
      : _ldt(ldt), _dt(dt), _ode(ode)
    {
      // Do nothing
    }

    // Constructor
    ODESolver (double ldt=0.0, double dt=0.0) : _ldt(ldt), _dt(dt), _ode(0)
    {
      // Do nothing
    }

    // Destructor
    virtual ~ODESolver () { /* Do nothing */ }

    // Attach ODE to solver
    virtual void attach(goss::ODE* ode) {_ode = ode;}

    // Step solver an interval of time forward
    virtual void forward(double* y, double t, double interval) = 0;

    // The size of the ODE
    inline uint ode_size() const { return _ode->size(); }

    // Return the ODE
    inline const goss::ODE* get_ode() const { return _ode; }

    // Return the internal time step
    inline double internal_time_step() const { return _ldt; }

    // Set the internal time step
    inline void set_internal_time_step(double ldt) { _ldt = ldt; }

  protected:
    
    // Local time step.
    double _ldt;

    // Variable local time step.
    double _dt;
    
    // Pointer to ode
    goss::ODE* _ode;

  };

}
#endif
