// Copyright (C) 2006-2012 Ola Skavhaug
//
// This file is part of GOSS.
//
// GOSS is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GOSS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GOSS. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Johan Hake 2012

#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdio>

#include "GRL2.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL2::GRL2() : ODESolver(0.0, 0.0), y0(0), a(0), b(0), linear_terms(0), 
	       delta(1.0e-8), nbytes(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
GRL2::GRL2(boost::shared_ptr<ODE> ode) : ODESolver(0.0, 0.0),
	                y0(0), a(0), b(0), linear_terms(0), 
	                delta(1.0e-8), nbytes(0)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
GRL2::GRL2(const GRL2& solver) : ODESolver(solver),
                                 y0(solver.y0), 
                                 a(solver.a), 
                                 b(solver.b), 
                                 linear_terms(solver.linear_terms),
                                 delta(solver.delta)
{
  assert(_ode);
  
  printf("GRL2 constructor!\n");
  // Get what terms are linear
  _ode->linear_terms(&linear_terms[0]);
}
//-----------------------------------------------------------------------------
GRL2::~GRL2()
{ 
  // Do nothing
}
//-----------------------------------------------------------------------------
void GRL2::attach(boost::shared_ptr<ODE> ode)
{
  // Attach ode using base class
  ODESolver::attach(ode);

  // Initalize memory
  y0.resize(num_states(), 0.0);
  a.resize(num_states(), 0.0);
  b.resize(num_states(), 0.0);
  linear_terms.resize(num_states());
  
  // Get what terms are linear
  _ode->linear_terms(&linear_terms[0]);
  
  nbytes = num_states()*sizeof(double);

}
//-----------------------------------------------------------------------------
void GRL2::forward(double* y, double t, double interval)
{
  
  assert(_ode);
  
  // Local timestep
  const double dt = interval;

  // Copy start conditions
  for (uint i = 0; i < num_states(); ++i) 
    y0[i] = y[i];
  //std::memcpy(&y0[0], y, nbytes); 

  // First step

  // Evaluate full right hand side
  _ode->eval(y, t, &a[0]);

  // Exact derivatives for linear terms
  _ode->linear_derivatives(y, t, &b[0]);  
  
  for (uint i = 0; i < num_states(); ++i) 
  { 
    // Numerical differentiation
    if (linear_terms[i] == 0) 
    {
      y[i] += delta; 
      b[i] = (_ode->eval(i, y, t) - a[i])/delta;  // Component i derivative
      y[i] -= delta;                       // Restore state i
    }
  }

  for (uint i = 0; i < num_states(); ++i) 
    y[i] += (std::fabs(b[i]) > delta) ? a[i]/b[i]*(std::exp(b[i]*dt*0.5) - 1.0) : 
      a[i]*dt*0.5;

  // Second step
  
  // Exact derivatives for linear terms
  //_ode->linear_derivatives(y, t, b);

  // Local variable to store comp i
  double yi;       
  for (uint i = 0; i < num_states(); ++i) 
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

  for (uint i = 0; i < num_states(); ++i) 
    y[i] = (std::fabs(b[i]) > delta) ? y0[i] + a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) :
      y0[i] + a[i]*dt;
}
//-----------------------------------------------------------------------------
