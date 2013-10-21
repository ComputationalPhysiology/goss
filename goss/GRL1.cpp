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
#include <cmath>
#include <cstdlib>

#include "log.h"
#include "GRL1.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL1::GRL1() : ODESolver(0.0, 0.0), a(0), b(0), linear_terms(0), 
	       delta(1.0e-8)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
GRL1::GRL1(boost::shared_ptr<ODE> ode) : ODESolver(0.0, 0.0), a(0), b(0), 
                                         linear_terms(0), delta(1.0e-8)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
GRL1::GRL1(const GRL1& solver) : ODESolver(solver), 
                                 a(solver.num_states()), 
                                 b(solver.num_states(), 0.0), 
                                 linear_terms(solver.num_states()),
				 delta(solver.delta)
{
  if (_ode)
    // Get what terms are linear
    _ode->linear_terms(&linear_terms[0]);
}
//-----------------------------------------------------------------------------
GRL1::~GRL1()
{ 
  // Do nothing
}

//-----------------------------------------------------------------------------
void GRL1::attach(boost::shared_ptr<ODE> ode)
{
  // Attach ode using base class
  ODESolver::attach(ode);
  
  if (ode->is_dae())
    goss_error("GRL1.cpp",
	       "attaching ode",
	       "cannot integrate a DAE ode with an explicit solver.");

  // Initalize memory
  a.resize(num_states());
  b.resize(num_states(), 0.0);
  linear_terms.resize(num_states());
  
  // Get what terms are linear
  _ode->linear_terms(&linear_terms[0]);
}

//-----------------------------------------------------------------------------
void GRL1::forward(double* y, double t, double interval)
{

  assert(_ode);

  // Local timestep
  const double dt = interval;

  // Evaluate full right hand side
  _ode->eval(y, t, &a[0]);              

  // Exact derivatives for linear terms 
  _ode->linear_derivatives(y, t, &b[0]);

  for (uint i = 0; i < num_states(); ++i) 
  { 
    // Numerical differentiation for non linear terms
    if (linear_terms[i] == 0) 
    {      
      y[i] += delta; 
      
      // Component i derivative
      b[i] = (_ode->eval(i, y, t) - a[i])/delta;  
      y[i] -= delta;				        // Restore state i
    }
  }

  // Integrate linear terms exactly
  for (uint i = 0; i < num_states(); ++i) 
    y[i] += (std::fabs(b[i]) > delta) ? a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) : a[i]*dt;

}
//-----------------------------------------------------------------------------
