#include <cmath>

#include "RK4.h"

using namespace goss;

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
RK4::~RK4 ()
{
  if (k1) delete[] k1;
  if (k2) delete[] k2;
  if (k3) delete[] k3;
  if (k4) delete[] k4;
  if (tmp) delete[] tmp;
}
//-----------------------------------------------------------------------------
void RK4::attach(ODE* ode)
{
  _ode = ode;
  if (k1) delete[] k1;
  if (k2) delete[] k2;
  if (k3) delete[] k3;
  if (k4) delete[] k4;
  if (tmp) delete[] tmp;

  k1  = new double[ode_size()];
  k2  = new double[ode_size()];
  k3  = new double[ode_size()];
  k4  = new double[ode_size()];
  tmp = new double[ode_size()];
}
//-----------------------------------------------------------------------------
void RK4::forward(double* y, double t, double interval) 
{
  // Calculate number of steps and size of timestep based on _ldt
  const ulong nsteps = _ldt > 0 ? std::ceil(interval/_ldt - 1.0E-12) : 1;
  const double dt = interval/nsteps;

  // Local time
  double lt = t;
  for (ulong j = 0; j < nsteps; ++j) 
  {
    // Evaluate rhs and calculate intermediate derivatives
    _ode->eval(y    , lt     , k1);

    // Explicit Euler step
    axpy(tmp, y, 0.5*dt, k1);

    _ode->eval(tmp, lt+0.5*dt, k2);
    axpy(tmp, y, 0.5*dt, k2);

    _ode->eval(tmp, lt+0.5*dt, k3);
    axpy(tmp, y, dt, k3);
    
    _ode->eval(tmp, lt+    dt, k4);

    for (uint i = 0; i < ode_size(); ++i)
      y[i] += dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
    
    // Update local time
    lt += dt;
  }
}
//-----------------------------------------------------------------------------
