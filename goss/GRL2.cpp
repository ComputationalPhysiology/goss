// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-19

#include <cassert>
#include <cstring>
#include <cmath>

#include "GRL2.h"
#include "LinearizedODE.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL2::GRL2() : ODESolver(0.0, 0.0), _lode(0), y0(0), a(0), b(0), linear_terms(0), 
	       delta(1.0e-8), nbytes(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
GRL2::GRL2(ODE* ode) : ODESolver(0.0, 0.0), _lode(0), y0(0), a(0), b(0), 
		       linear_terms(0), delta(1.0e-8), nbytes(0)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
GRL2::GRL2(const GRL2& solver) : ODESolver(solver), _lode(0), 
				 y0(new double[solver.ode_size()]), 
				 a(new double[solver.ode_size()]), 
				 b(new double[solver.ode_size()]), 
				 linear_terms(new uint[solver.ode_size()]),
				 delta(solver.delta)
{
  // Store Linearized ODE
  _lode = dynamic_cast<LinearizedODE*>(_ode.get());
  assert(_lode);

  // Get what terms are linear
  _lode->linear_terms(linear_terms.get());
}
//-----------------------------------------------------------------------------
GRL2::~GRL2()
{ 
  // Do nothing
}
//-----------------------------------------------------------------------------
void GRL2::attach(ODE* ode)
{
  // Attach ode using base class
  ODESolver::attach(ode);

   // Store Linearized ODE
  _lode = dynamic_cast<LinearizedODE*>(ode);
  assert(_lode);
  
  // Initalize memory
  y0.reset(new double[ode_size()]);
  a.reset(new double[ode_size()]);
  b.reset(new double[ode_size()]);
  linear_terms.reset(new uint[ode_size()]);
  std::fill(b.get(), b.get()+ode_size(), static_cast<double>(0));
  
  // Get what terms are linear
  _lode->linear_terms(linear_terms.get());
  
  nbytes  = ode_size()*sizeof(double);
  std::fill(a.get(), a.get() + ode_size(), 0.0);
  std::fill(b.get(), b.get() + ode_size(), 0.0);
  std::fill(y0.get(), y0.get() + ode_size(), 0.0);
}
//-----------------------------------------------------------------------------
void GRL2::forward(double* y, double t, double interval)
{
  
  assert(_lode);
  
  // Local timestep
  const double dt = interval;

  // Copy start conditions
  std::memcpy(y0.get(), y, nbytes); 

  // First step

  // Evaluate full right hand side
  _lode->eval(y, t, a.get());

  // Exact derivatives for linear terms
  _lode->linear_derivatives(y, t, b.get());  
  
  for (uint i = 0; i < ode_size(); ++i) 
  { 
    // Numerical differentiation
    if (linear_terms[i] == 0) 
    {
      y[i] += delta; 
      b[i] = (_ode->eval(i, y, t) - a[i])/delta;  // Component i derivative
      y[i] -= delta;				        // Restore state i
    }
  }

  for (uint i = 0; i < ode_size(); ++i) 
    y[i] += (std::fabs(b[i]) > delta) ? a[i]/b[i]*(std::exp(b[i]*dt*0.5) - 1.0) : 
      a[i]*dt*0.5;

  // Second step
  
  // Exact derivatives for linear terms
  //_lode->linear_derivatives(y, t, b);

  // Local variable to store comp i
  double yi;					        
  for (uint i = 0; i < ode_size(); ++i) 
  {        
    // Store original value of comp i
    yi = y[i];

    // Use y0[i] at diagonal
    y[i] = y0[i];

    // Evaluate right hand side after changing comp i    
    a[i] = _ode->eval(i, y, t);

//    if (not linear_terms[i]) 
//    {  
      // Numerical differentiation for non-linear terms
      y[i] += delta; 

      // Component i derivative
      b[i] = (_ode->eval(i, y, t) - a[i])/delta;
//    }

    // Restore state i
    y[i] = yi;
  }

  for (uint i = 0; i < ode_size(); ++i) 
    y[i] = (std::fabs(b[i]) > delta) ? y0[i] + a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) :
      y0[i] + a[i]*dt;
}
//-----------------------------------------------------------------------------
