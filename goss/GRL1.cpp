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
#include "Timer.h"
#include "GRL1.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL1::GRL1() : ODESolver(), _a(0), _b(0), _delta(1.0e-8)
{
  parameters.rename("GRL1");
}
//-----------------------------------------------------------------------------
GRL1::GRL1(boost::shared_ptr<ODE> ode) : ODESolver(), _a(0), _b(0), 
                                         _delta(1.0e-8)
{
  parameters.rename("GRL1");
  attach(ode);
}
//-----------------------------------------------------------------------------
GRL1::GRL1(const GRL1& solver) : ODESolver(solver), 
                                 _a(solver.num_states()), 
                                 _b(solver.num_states(), 0.0), 
				 _delta(solver._delta)
{
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
  _a.resize(num_states());
  _b.resize(num_states(), 0.0);
}

//-----------------------------------------------------------------------------
void GRL1::forward(double* y, double t, double dt)
{

  Timer _timer("GRL1 one step");

  assert(_ode);

  // Calculate number of steps and size of timestep based on _ldt
  const double ldt_0 = parameters["ldt"];
  const ulong nsteps = ldt_0 > 0 ? std::ceil(dt/ldt_0 - 1.0E-12) : 1;
  const double ldt = dt/nsteps;
  
  // Local time
  double lt = t;

  for (ulong step = 0; step < nsteps; ++step)
  {

    // Evaluate full right hand side
    _ode->linearized_eval(y, lt, _b.data(), _a.data());
    //_ode->eval(y, lt, a.data());              
  
    // Exact derivatives for linear terms 
    //_ode->linear_derivatives(y, lt, b.data());
    //
    //for (uint i = 0; i < num_states(); ++i) 
    //{ 
    //  // Numerical differentiation for non linear terms
    //  if (!_ode->linear_term(i))
    //  {      
    //    y[i] += delta; 
    //    
    //    // Component i derivative
    //    b[i] = (_ode->eval(i, y, lt) - a[i])/delta;
    //    y[i] -= delta;				   // Restore state i
    //  }
    //}
  
    // Integrate linear terms exactly
    for (uint i = 0; i < num_states(); ++i) 
      y[i] += (std::fabs(_b[i]) > _delta) ? _a[i]/_b[i]*(std::exp(_b[i]*ldt) - 1.0) : _a[i]*ldt;
  
      // Increase time
    lt += ldt;
  }
}
//-----------------------------------------------------------------------------
