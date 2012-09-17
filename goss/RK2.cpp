#include <cstdlib>
#include "RK2.h"

using namespace goss;

//-----------------------------------------------------------------------------
RK2::RK2() : ODESolver(), k1(0), tmp(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RK2::RK2(ODE* ode) : ODESolver(), k1(0), tmp(0)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
RK2::~RK2 ()
{
  if (k1) delete[] k1;
  if (tmp) delete[] tmp;
}
//-----------------------------------------------------------------------------
void RK2::attach(ODE* ode)
{
  this->_ode = ode;
  if (k1) delete[] k1;
  if (tmp) delete[] tmp;

  k1  = new double[ode_size()];
  tmp = new double[ode_size()];
}
//-----------------------------------------------------------------------------
void RK2::forward(double* y, double t, double dt) 
{
  // Initial eval
  _ode->eval(y, t, k1);

  // Explicit Euler step to find the midpoint solution
  axpy(tmp, y, 0.5*dt, k1);
  
  // Evaluate derivative at midpoint
  _ode->eval(tmp, t+0.5*dt, k1);

  // Use midpoint derivative for explicit Euler step
  for (uint i = 0; i < ode_size(); ++i)
    y[i] += dt*k1[i];

}
//-----------------------------------------------------------------------------
