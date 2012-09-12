#include "ExplicitEuler.h"

using namespace gossol;

//-----------------------------------------------------------------------------
ExplicitEuler:: ExplicitEuler(ODE* ode, double ldt)
//-----------------------------------------------------------------------------
  : ODESolver(ode, ldt), dFdt(0)
{ 
  attach(ode);
} 

//-----------------------------------------------------------------------------
ExplicitEuler:: ExplicitEuler(double ldt) 
//-----------------------------------------------------------------------------
  : ODESolver(ldt), dFdt(0)
{
  // Do nothing
}
     

//-----------------------------------------------------------------------------
void  ExplicitEuler:: attach(ODE* ode)
//-----------------------------------------------------------------------------
{
  free(dFdt);
  this->ode = ode;
  dFdt = static_cast<double*>(malloc(ode->size()*sizeof*dFdt));
}


//-----------------------------------------------------------------------------
void ExplicitEuler:: forward(double* y, double t, double interval) {
//-----------------------------------------------------------------------------
  int n_steps = 1;
 
  if (ldt>0) {
    // Calculates a local time step that is <= ldt,
    // and that divides dt into equally sized steps:
    n_steps = (int) ceil(interval/ldt - 1.0E-12); 
    dt = interval/((double)n_steps);
  } else {
    // No internal time step chosen, use interval instead
    dt = interval;
  }
 
  double lt = t;
  for (int j = 0; j<n_steps; ++j){
    ode->eval(y, lt, dFdt);
    for (int i=0; i < ode->size(); ++i)
      y[i] += dt*dFdt[i];
    lt += dt;
  }
}


