// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2007-07-09

#include <iostream>
#include <cstdlib>
#include "GRL1.h"

using namespace gossol;

//-----------------------------------------------------------------------------
GRL1::GRL1()
//-----------------------------------------------------------------------------
  : ODESolver(0.0, 0.0), a(0), b(0), linear_terms(0), delta(1.0e-8)
{
  // Do nothing
}

//-----------------------------------------------------------------------------
GRL1:: GRL1(gossol::ODE* ode)
//-----------------------------------------------------------------------------
  : ODESolver(ode, 0.0), a(0), b(0), linear_terms(0), delta(1.0e-8)
{
  attach(ode);
}

//-----------------------------------------------------------------------------
GRL1::~GRL1()
//-----------------------------------------------------------------------------
{ 
  free(a); free(b); 
  free(linear_terms);
}

//-----------------------------------------------------------------------------
void GRL1::attach(ODE* ode)
//-----------------------------------------------------------------------------
{
  free(a); free(b);
  this->ode = ode;
  n = ode->size();
  a = static_cast<double*>(malloc(n*sizeof*a));
  b = static_cast<double*>(malloc(n*sizeof*b));
  std::fill(a, a+n, static_cast<double>(0));
  std::fill(b, b+n, static_cast<double>(0));
  linear_terms = static_cast<int*>(malloc(n*sizeof*linear_terms));
  ode->linearTerms(linear_terms);
}

//-----------------------------------------------------------------------------
void GRL1::forward(double* y, double t, double dt)
//-----------------------------------------------------------------------------
{

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
    y[i] += (fabs(b[i]) > delta) ? a[i]/b[i]*(exp(b[i]*dt) - 1.0) : a[i]*dt;
}
