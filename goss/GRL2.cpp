// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2007-10-29

#include <iostream>
#include <cstdlib>
#include <cstring>
#include "GRL2.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL2::GRL2()
//-----------------------------------------------------------------------------
  : ODESolver(0.0, 0.0), y0(0), a(0), b(0), linear_terms(0), delta(1.0e-8)
{
  // Do nothing
}

//-----------------------------------------------------------------------------
GRL2:: GRL2(goss::ODE* ode)
//-----------------------------------------------------------------------------
  : ODESolver(ode, 0.0), y0(0), a(0), b(0), linear_terms(0), delta(1.0e-8)
{
  attach(ode);
}

//-----------------------------------------------------------------------------
GRL2::~GRL2()
//-----------------------------------------------------------------------------
{ 
  free(y0);
  free(a); 
  free(b); 
  free(linear_terms);
}

//-----------------------------------------------------------------------------
void GRL2::attach(ODE* ode)
//-----------------------------------------------------------------------------
{
  free(a); free(b);
  this->ode = ode;
  n = ode->size();
  y0 = static_cast<double*>(malloc(n*sizeof*y0));
  a =  static_cast<double*>(malloc(n*sizeof*a));
  b =  static_cast<double*>(malloc(n*sizeof*b));
  linear_terms = static_cast<int*>(malloc(n*sizeof*linear_terms));
  ode->linearTerms(linear_terms);
  std::fill(a, a+n, static_cast<double>(0));
  std::fill(b, b+n, static_cast<double>(0));
  std::fill(y0, y0+n, static_cast<double>(0));
}

//-----------------------------------------------------------------------------
void GRL2::forward(double* y, double t, double dt)
//-----------------------------------------------------------------------------
{
  double yi;					        // Local variable to store comp i
  memcpy(y0, y, n*sizeof*y0); // Copy start conditions

  // First step
  ode->eval(y, t, a);               // Evaluate full right hand side
  ode->linearDerivatives(y, t, b);  // Exact derivatives for linear terms
  for (int i=0; i<n; ++i) { 
    if (linear_terms[i] == 0) {      // Numerical differentiation
      y[i] += delta; 
      b[i] = (ode->eval(i, y, t) - a[i])/delta;  // Component i derivative
      y[i] -= delta;				        // Restore state i
    }
  }
  for (int i=0; i<n; ++i) 
    y[i] += (fabs(b[i]) > delta) ? a[i]/b[i]*(exp(b[i]*dt*0.5) - 1.0) : a[i]*dt*0.5;

  // Second step
  //ode->linearDerivatives(y, t, b);  // Exact derivatives for linear terms
  for (int i=0; i<n; ++i) {        
    yi = y[i];				              // Store original value of comp i
    y[i] = y0[i];                   // Use y0[i] at diagonal
    a[i] = ode->eval(i, y, t);      // Evaluate right hand side after changing comp i 

//    if (not linear_terms[i]) {  // Numerical differentiation for non-linear terms
      y[i] += delta; 
      b[i] = (ode->eval(i, y, t) - a[i])/delta;  // Component i derivative
//    }
    y[i] = yi;                  // Restore state i
  }
  for (int i=0; i<n; ++i)
    y[i] = (fabs(b[i]) > delta) ? y0[i] + a[i]/b[i]*(exp(b[i]*dt) - 1.0) : y0[i] + a[i]*dt;
}

