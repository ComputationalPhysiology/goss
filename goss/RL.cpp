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

#include "RL.h"

using namespace goss;

//-----------------------------------------------------------------------------
RL::RL() : ODESolver(0.0, 0.0), a(0), b(0), linear_terms(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
RL::RL(boost::shared_ptr<ODE> ode) : ODESolver(0.0, 0.0), a(0), b(0), 
				     linear_terms(0)
{
  attach(ode);
}
//-----------------------------------------------------------------------------
RL::RL(const RL& solver) : ODESolver(solver), 
			   a(solver.num_states(), 0.0), 
			   b(solver.num_states(), 0.0), 
			   linear_terms(solver.num_states())
{
  if (_ode)
    
    // Get what terms are linear
    _ode->linear_terms(&linear_terms[0]);
}
//-----------------------------------------------------------------------------
RL::~RL()
{ 
  // Do nothing
}
//-----------------------------------------------------------------------------
void RL::attach(boost::shared_ptr<ODE> ode)
{
  // Attach ode using base class
  ODESolver::attach(ode);
  
  // Initalize memory
  a.resize(num_states(), 0.0);
  b.resize(num_states(), 0.0);
  linear_terms.resize(num_states());
  
  // Get what terms are linear
  _ode->linear_terms(&linear_terms[0]);
}
//-----------------------------------------------------------------------------
void RL::forward(double* y, double t, double interval)
{

  assert(_ode);

  if (_ode->is_dae())
    goss_error("RL.cpp",
	       "forwarding ode",
	       "cannot integrate a DAE ode with an explicit solver.");

  // Local timestep
  const double dt = interval;

  // Evaluate full right hand side
  _ode->eval(y, t, &a[0]);

  // Exact derivatives for linear terms
  _ode->linear_derivatives(y, t, &b[0]);

  // Integrate linear terms exactly
  for (uint i = 0; i < num_states(); ++i) 
    y[i] += (linear_terms[i]==1) ? a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) : a[i]*dt;

}
//-----------------------------------------------------------------------------
