// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-19

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "GRL1.h"
#include "LinearizedODE.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL1::GRL1() : ODESolver(0.0, 0.0), _lode(0), a(0), b(0), linear_terms(0), 
	       delta(1.0e-8)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
GRL1::GRL1(ODE* ode) : ODESolver(0.0, 0.0), _lode(0), a(0), b(0), 
		       linear_terms(0), delta(1.0e-8)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
GRL1::GRL1(const GRL1& solver) : ODESolver(solver), _lode(0), 
				 a(new double[solver.num_states()]), 
				 b(new double[solver.num_states()]), 
				 linear_terms(new uint[solver.num_states()]),
				 delta(solver.delta)
{
  // Store Linearized ODE
  _lode = dynamic_cast<LinearizedODE*>(_ode.get());
  assert(_lode);

  // Get what terms are linear
  _lode->linear_terms(linear_terms.get());
}
//-----------------------------------------------------------------------------
GRL1::~GRL1()
{ 
  // Do nothing
}

//-----------------------------------------------------------------------------
void GRL1::attach(ODE* ode)
{
  // Attach ode using base class
  ODESolver::attach(ode);
  
  // Store Linearized ODE
  _lode = dynamic_cast<LinearizedODE*>(ode);
  assert(_lode);
  
  // Initalize memory
  a.reset(new double[num_states()]);
  b.reset(new double[num_states()]);
  linear_terms.reset(new uint[num_states()]);
  std::fill(b.get(), b.get()+num_states(), static_cast<double>(0));
  
  // Get what terms are linear
  _lode->linear_terms(linear_terms.get());
}

//-----------------------------------------------------------------------------
void GRL1::forward(double* y, double t, double interval)
{

  assert(_lode);

  // Local timestep
  const double dt = interval;

  // Evaluate full right hand side
  _lode->eval(y, t, a.get());              

  // Exact derivatives for linear terms 
  _lode->linear_derivatives(y, t, b.get());

  for (uint i = 0; i < num_states(); ++i) 
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
  for (uint i = 0; i < num_states(); ++i) 
    y[i] += (std::fabs(b[i]) > delta) ? a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) : a[i]*dt;

}
