// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-17

#include <cmath>
#include <cstdlib>
#include "GRL1.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL1::GRL1() : ODESolver(0.0, 0.0), _lode(0), a(0), b(0), linear_terms(0), 
	       delta(1.0e-8)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
GRL1::GRL1(LinearizedODE* ode) : ODESolver(ode, 0.0), _lode(0), a(0), b(0), 
				 linear_terms(0), delta(1.0e-8)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
GRL1::~GRL1()
{ 
  if (a) delete[] a;
  if (b) delete[] b;
  if (linear_terms) delete[] linear_terms;
}

//-----------------------------------------------------------------------------
void GRL1::attach(LinearizedODE* ode)
{
  if (a) delete[] a;
  if (b) delete[] b;
  if (linear_terms) delete[] linear_terms;

  // Store Linearized and ordinary ODE
  _ode = ode;
  _lode = ode;
  
  // Initalize memory
  a = new double[ode_size()];
  b = new double[ode_size()];
  linear_terms = new uint[ode_size()];
  std::fill(b, b+ode_size(), static_cast<double>(0));
  
  // Get what terms are linear
  _lode->linear_terms(linear_terms);
}

//-----------------------------------------------------------------------------
void GRL1::forward(double* y, double t, double interval)
//-----------------------------------------------------------------------------
{

  // Local timestep
  const double dt = interval;

  // Evaluate full right hand side
  _lode->eval(y, t, a);              

  // Exact derivatives for linear terms 
  _lode->linear_derivatives(y, t, b);

  for (uint i = 0; i < ode_size(); ++i) 
  { 
    // Numerical differentiation for non linear terms
    if (linear_terms[i] == 0) 
    {      
      y[i] += delta; 
      
      // Component i derivative
      b[i] = (_lode->eval(i, y, t) - a[i])/delta;  
      y[i] -= delta;				        // Restore state i
    }
  }

  // Integrate linear terms exactly
  for (uint i = 0; i < ode_size(); ++i) 
    y[i] += (std::fabs(b[i]) > delta) ? a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) : a[i]*dt;

}
