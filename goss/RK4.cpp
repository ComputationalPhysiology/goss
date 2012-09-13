#include <cstdlib>
#include "RK4.h"

using namespace goss;

//-----------------------------------------------------------------------------
RK4::RK4( double _ldt)
//-----------------------------------------------------------------------------
{
  ldt = _ldt;
}

//-----------------------------------------------------------------------------
RK4::RK4(ODE *ode_, double _ldt)
//-----------------------------------------------------------------------------
{
  attach(ode_);
  ldt = _ldt; 
}

//-----------------------------------------------------------------------------
RK4::~RK4 ()
//-----------------------------------------------------------------------------
{
  free(k1); free(k2);
  free(k3); free(k4);
  free(tmp);
}

//-----------------------------------------------------------------------------
void  RK4:: attach(ODE* ode_)
//-----------------------------------------------------------------------------
{
  ode = ode_;
  k1  = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  k2  = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  k3  = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  k4  = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  tmp = static_cast<double*>(malloc(sizeof(double)*ode->size()));
}


//-----------------------------------------------------------------------------
void RK4:: forward(double* y, double t, double interval) 
//-----------------------------------------------------------------------------
{
  int n_steps = 1;

  if ( ldt>0 ) {
    // Calculates a local time step that is <= ldt,
    // and that divides dt into equally sized steps:

    n_steps = (int) ceil(interval/ldt -1E-12); 
    dt =interval/((double)n_steps); 
  }

  double lt = t;
  for (int j = 0; j<n_steps; ++j) {
    ode->eval(y    , lt     , k1);
    add(tmp, y, 0.5*dt, k1);
    ode->eval(tmp, lt+0.5*dt, k2);
    add(tmp, y, 0.5*dt, k2);
    ode->eval(tmp, lt+0.5*dt, k3);
    add(tmp, y, dt, k3);
    ode->eval(tmp, lt+    dt, k4);

    for (int i=0; i < ode->size(); ++i)
      y[i] += dt*( k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i] )/6.0;
    lt += dt;
  }
}


