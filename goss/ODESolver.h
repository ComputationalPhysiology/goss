#ifndef ODESolver_h_IS_INCLUDED
#define ODESolver_h_IS_INCLUDED

#include <boost/shared_ptr.hpp>
#include <iostream>

#include "ODE.h"
#include "types.h"

namespace goss 
{

  class ODESolver {

  public:

    // FIXME: Remove dt as variable...
    // Default Constructor
    ODESolver () : _ldt(-1.0), _dt(0.0), _ode(static_cast<ODE*>(0))
    {
      // Do nothing
    }

    // Constructor
    ODESolver (double ldt, double dt=0.0) : _ldt(ldt), _dt(dt), 
					    _ode(static_cast<ODE*>(0))
    {
      // Do nothing
    }

    // Copy constructor (Uses default copy constructor of ODE)
    ODESolver (const ODESolver& solver) : _ldt(solver._ldt), _dt(solver._dt), 
					  _ode(solver._ode->copy())
    {
      // Do nothing
    }

    // Destructor
    virtual ~ODESolver () 
    { 
      // Do nothing
    }

    // Return a copy of itself
    virtual boost::shared_ptr<ODESolver> copy() const = 0;

    // Attach ODE and reset solver 
    virtual void attach(boost::shared_ptr<ODE> ode) 
    { _ode = ode; reset();}

    // Reset solver 
    virtual void reset() { /* Do nothing */ }

    // Step solver an interval of time forward
    virtual void forward(double* y, double t, double interval) = 0;

    // The size of the ODE
    inline uint num_states() const { return _ode->num_states(); }

    // Return the ODE (const version)
    inline const boost::shared_ptr<ODE> get_ode() const { return _ode; }

    // Return the ODE
    inline boost::shared_ptr<ODE> get_ode() { return _ode; }

    // Return the internal time step
    inline double get_internal_time_step() const { return _ldt; }

    // Set the internal time step
    inline void set_internal_time_step(double ldt) { _ldt = ldt; }

    // Return true if the Solver is adaptive
    virtual bool is_adaptive() const { return false; }

  protected:
    
    // Local time step.
    double _ldt;

    // Variable local time step.
    double _dt;
    
    // Shared pointer to ode 
    boost::shared_ptr<ODE> _ode;

  };

}
#endif
