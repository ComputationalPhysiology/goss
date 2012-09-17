#include <cassert>

#include "ExplicitEuler.h"

using namespace goss;

//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(ODE* ode, double ldt) : ODESolver(ode, ldt), _dFdt(0)
{ 
  attach(ode);
} 
//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(double ldt) : ODESolver(ldt), _dFdt(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void  ExplicitEuler::attach(ODE* ode)
{
  // Attach ODE
  _ode = ode;

  // Clean
  if (_dFdt)
    delete[] _dFdt;

  // Create memory for derivative evaluation
  _dFdt = new double[ode->size()];

}
//-----------------------------------------------------------------------------
void ExplicitEuler::forward(double* y, double t, double interval) 
{

  assert(_ode);
  
  // Calculate number of steps and size of timestep based on _ldt
  const ulong nsteps = _ldt > 0 ? ceil(interval/_ldt - 1.0E-12) : 1;
  const double dt = interval/nsteps;

  // Local time
  double lt = t;
  for (ulong step = 0; step < nsteps; ++step)
  {
    // Evaluate rhs
    _ode->eval(y, lt, _dFdt);

    // Update states
    for (uint i = 0; i < _ode->size(); ++i)
      y[i] += dt*_dFdt[i];

    // Increase time
    lt += dt;

  }

}
