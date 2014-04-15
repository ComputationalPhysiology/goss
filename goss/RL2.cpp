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

#include "log.h"
#include "Timer.h"
#include "RL2.h"

using namespace goss;

//-----------------------------------------------------------------------------
RL2::RL2() : ODESolver(), _y2(0), _a(0), _b(0)
{
  parameters.rename("RL2");
}
//-----------------------------------------------------------------------------
RL2::RL2(boost::shared_ptr<ODE> ode) : ODESolver(),
  _y2(0), _a(0), _b(0)
{
  parameters.rename("RL2");
  attach(ode);
}
//-----------------------------------------------------------------------------
RL2::RL2(const RL2& solver) : ODESolver(solver),
                                 _y2(solver._y2), 
                                 _a(solver._a), 
                                 _b(solver._b) 
{
}
//-----------------------------------------------------------------------------
RL2::~RL2()
{ 
  // Do nothing
}
//-----------------------------------------------------------------------------
void RL2::attach(boost::shared_ptr<ODE> ode)
{
  // Attach ode using base class
  ODESolver::attach(ode);

  if (ode->is_dae())
    goss_error("RL2.cpp",
	       "attaching ode",
	       "cannot integrate a DAE ode with an explicit solver.");

  // Initalize memory
  _y2.resize(num_states(), 0.0);
  _a.resize(num_states(), 0.0);
  _b.resize(num_states(), 0.0);
}
//-----------------------------------------------------------------------------
void RL2::one_step(double* y2, double* y, double* y0, double t, double dt)
{
  Timer _timer("RL2 one step");

  // Evaluate full right hand side
  _ode->linearized_eval(y, t, _b.data(), _a.data());
  //_ode->eval(y, t, _a.data());
  //
  //// Exact derivatives for linear terms
  //_ode->linear_derivatives(y, t, _b.data());
  
  for (uint i = 0; i < num_states(); ++i) 
    y2[i] = _ode->linear_term(i) ? y0[i] + _a[i]/_b[i]*(std::exp(_b[i]*dt) - 1.0) :
      y0[i] + _a[i]*dt;
}
//-----------------------------------------------------------------------------
void RL2::forward(double* y, double t, double dt)
{
  
  assert(_ode);
  
  // First step
  one_step(_y2.data(), y, y, t, dt*0.5);

  one_step(y, _y2.data(), y, t, dt);

}
//-----------------------------------------------------------------------------
