// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-19

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "RL.h"
#include "LinearizedODE.h"

using namespace goss;

//-----------------------------------------------------------------------------
RL::RL() : ODESolver(0.0, 0.0), _lode(0), a(0), b(0), linear_terms(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RL::RL(ODE* ode) : ODESolver(0.0, 0.0), _lode(0), a(0), b(0), 
		   linear_terms(0)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
RL::RL(const RL& solver) : ODESolver(solver), _lode(0), 
			   a(new double[solver.num_states()]), 
			   b(new double[solver.num_states()]), 
			   linear_terms(new uint[solver.num_states()])
{
  // Store Linearized ODE
  _lode = dynamic_cast<LinearizedODE*>(_ode.get());
  assert(_lode);

  // Get what terms are linear
  _lode->linear_terms(linear_terms.get());
}
//-----------------------------------------------------------------------------
RL::~RL()
{ 
  // Do nothing
}
//-----------------------------------------------------------------------------
void RL::attach(ODE* ode)
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
void RL::forward(double* y, double t, double interval)
{

  assert(_lode);

  // Local timestep
  const double dt = interval;

  // Evaluate full right hand side
  _lode->eval(y, t, a.get());

  // Exact derivatives for linear terms
  _lode->linear_derivatives(y, t, b.get());

  // Integrate linear terms exactly
  for (uint i = 0; i < num_states(); ++i) 
    y[i] += (linear_terms[i]==1) ? a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) : a[i]*dt;

}
//-----------------------------------------------------------------------------
