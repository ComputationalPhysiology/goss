
#include <cassert>
#include <cmath>

#include "RK2.h"

using namespace goss;

//-----------------------------------------------------------------------------
RK2::RK2() : ODESolver(), k1(0), tmp(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RK2::RK2(double ldt) : ODESolver(ldt), k1(0), tmp(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RK2::RK2(ODE* ode, double ldt) : ODESolver(ldt), k1(0), tmp(0)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
RK2::RK2(const RK2& solver) : ODESolver(solver), k1(new double[solver.num_states()]),
			      tmp(new double[solver.num_states()])
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RK2::~RK2 ()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void RK2::attach(ODE* ode)
{
  
  ODESolver::attach(ode);

  k1.reset(new double[num_states()]);
  tmp.reset(new double[num_states()]);
}
//-----------------------------------------------------------------------------
void RK2::forward(double* y, double t, double interval) 
{

  assert(_ode);

  // Calculate number of steps and size of timestep based on _ldt
  const ulong nsteps = _ldt > 0 ? std::ceil(interval/_ldt - 1.0E-12) : 1;
  const double dt = interval/nsteps;

  // Local time
  double lt = t;
  for (ulong j = 0; j < nsteps; ++j) 
  {
    // Initial eval
    _ode->eval(y, lt, k1.get());

    // Explicit Euler step to find the midpoint solution
    axpy(tmp.get(), y, 0.5*dt, k1.get());
  
    // Evaluate derivative at midpoint
    _ode->eval(tmp.get(), lt+0.5*dt, k1.get());

    // Use midpoint derivative for explicit Euler step
    for (uint i = 0; i < num_states(); ++i)
      y[i] += dt*k1[i];

    // Update local time
    lt += dt;
  }
}
//-----------------------------------------------------------------------------
