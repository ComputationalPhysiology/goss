#include <cassert>
#include <cmath>

#include "RK4.h"

using namespace goss;

//-----------------------------------------------------------------------------
RK4::RK4() : ODESolver(), k1(0), k2(0), k3(0), k4(0), tmp(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RK4::RK4(double ldt) : ODESolver(ldt), k1(0), k2(0), k3(0), k4(0), tmp(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RK4::RK4(ODE *ode, double ldt) : ODESolver(ldt), k1(0), k2(0), k3(0), k4(0), tmp(0)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
RK4::RK4(const RK4& solver) : ODESolver(solver), k1(new double[solver.num_states()]),
			      k2(new double[solver.num_states()]),
			      k3(new double[solver.num_states()]),
			      k4(new double[solver.num_states()]),
			      tmp(new double[solver.num_states()])
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RK4::~RK4 ()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void RK4::attach(ODE* ode)
{

  // Attach ode using base class
  ODESolver::attach(ode);
  
  k1.reset(new double[num_states()]);
  k2.reset(new double[num_states()]);
  k3.reset(new double[num_states()]);
  k4.reset(new double[num_states()]);
  tmp.reset(new double[num_states()]);

}
//-----------------------------------------------------------------------------
void RK4::forward(double* y, double t, double interval) 
{

  assert(_ode);

  // Calculate number of steps and size of timestep based on _ldt
  const ulong nsteps = _ldt > 0 ? std::ceil(interval/_ldt - 1.0E-12) : 1;
  const double dt = interval/nsteps;

  // Local time
  double lt = t;
  for (ulong j = 0; j < nsteps; ++j) 
  {
    // Evaluate rhs and calculate intermediate derivatives
    _ode->eval(y, lt, k1.get());

    // Explicit Euler step
    axpy(tmp.get(), y, 0.5*dt, k1.get());

    _ode->eval(tmp.get(), lt + 0.5*dt, k2.get());
    axpy(tmp.get(), y, 0.5*dt, k2.get());

    _ode->eval(tmp.get(), lt + 0.5*dt, k3.get());
    axpy(tmp.get(), y, dt, k3.get());
    
    _ode->eval(tmp.get(), lt + dt, k4.get());

    for (uint i = 0; i < num_states(); ++i)
      y[i] += dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
    
    // Update local time
    lt += dt;
  }
}
//-----------------------------------------------------------------------------
