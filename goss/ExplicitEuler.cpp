#include <cassert>
#include <cmath>

#include "ExplicitEuler.h"

using namespace goss;

//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler() : ODESolver(), _dFdt(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(double ldt) : ODESolver(ldt), _dFdt(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(boost::shared_ptr<ODE> ode, double ldt) : 
  ODESolver(ldt), _dFdt(0)
{ 
  attach(ode);
} 
//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(const ExplicitEuler& solver) : 
  ODESolver(solver), _dFdt(new double[solver.num_states()])
{ 
  // Do nothing
} 
//-----------------------------------------------------------------------------
ExplicitEuler::~ExplicitEuler()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ExplicitEuler::attach(boost::shared_ptr<ODE> ode)
{
  // Attach ODE
  ODESolver::attach(ode);

  // Create memory for derivative evaluation
  _dFdt.reset(new double[num_states()]);

}
//-----------------------------------------------------------------------------
void ExplicitEuler::forward(double* y, double t, double interval) 
{

  assert(_ode);
  
  // Calculate number of steps and size of timestep based on _ldt
  const ulong nsteps = _ldt > 0 ? std::ceil(interval/_ldt - 1.0E-12) : 1;
  const double dt = interval/nsteps;

  // Local time
  double lt = t;
  for (ulong step = 0; step < nsteps; ++step)
  {
    // Evaluate rhs
    _ode->eval(y, lt, _dFdt.get());

    // Update states
    for (uint i = 0; i < num_states(); ++i)
      y[i] += dt*_dFdt[i];

    // Increase time
    lt += dt;

  }

}
