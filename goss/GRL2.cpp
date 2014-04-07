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
#include "GRL2.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL2::GRL2() : ODESolver(), _y2(0), _a(0), _b(0), _linear_terms(0), 
	       _delta(1.0e-8)
{
  parameters.rename("GRL2");
}
//-----------------------------------------------------------------------------
GRL2::GRL2(boost::shared_ptr<ODE> ode) : ODESolver(),
	                _y2(0), _a(0), _b(0), _linear_terms(0), 
	                _delta(1.0e-8)
{
  parameters.rename("GRL2");
  attach(ode);
}
//-----------------------------------------------------------------------------
GRL2::GRL2(const GRL2& solver) : ODESolver(solver),
                                 _y2(solver._y2), 
                                 _a(solver._a), 
                                 _b(solver._b), 
                                 _linear_terms(solver._linear_terms),
                                 _delta(solver._delta)
{
  
  if (_ode)
    // Get what terms are linear
    _ode->linear_terms(_linear_terms.data());
  
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

  if (ode->is_dae())
    goss_error("GRL2.cpp",
	       "attaching ode",
	       "cannot integrate a DAE ode with an explicit solver.");

  // Initalize memory
  _y2.resize(num_states(), 0.0);
  _a.resize(num_states(), 0.0);
  _b.resize(num_states(), 0.0);
  _linear_terms.resize(num_states());
  
  // Get what terms are linear
  _ode->linear_terms(_linear_terms.data());
  
}
//-----------------------------------------------------------------------------
void GRL2::one_step(double* y2, double* y, double* y0, double t, double dt)
{
  Timer _timer("GRL2 one step");

  // Evaluate full right hand side
  _ode->eval(y, t, _a.data());

  // Exact derivatives for linear terms
  _ode->linear_derivatives(y, t, _b.data());
  
  for (uint i = 0; i < num_states(); ++i) 
  { 
    // Numerical differentiation
    if (_linear_terms[i] == 0)
    {
      y[i] += _delta; 
      _b[i] = (_ode->eval(i, y, t) - _a[i])/_delta;  // Component i derivative
      y[i] -= _delta;                       // Restore state i
    }
  }

  for (uint i = 0; i < num_states(); ++i) 
    y2[i] = (std::fabs(_b[i]) > _delta) ? y0[i] + _a[i]/_b[i]*(std::exp(_b[i]*dt) - 1.0) :
      y0[i] + _a[i]*dt;
}
//-----------------------------------------------------------------------------
void GRL2::forward(double* y, double t, double dt)
{
  
  assert(_ode);
  
  // First step
  one_step(_y2.data(), y, y, t, dt*0.5);

  one_step(y, _y2.data(), y, t, dt);

  //// Evaluate full right hand side
  //_ode->eval(y, t, a.data());
  //
  //// Exact derivatives for linear terms
  //_ode->linear_derivatives(y, t, b.data());
  //
  //for (uint i = 0; i < num_states(); ++i) 
  //{ 
  //  // Numerical differentiation
  //  if (linear_terms[i] == 0) 
  //  {
  //    y[i] += delta; 
  //    b[i] = (_ode->eval(i, y, t) - a[i])/delta;  // Component i derivative
  //    y[i] -= delta;                       // Restore state i
  //  }
  //}
  //
  //for (uint i = 0; i < num_states(); ++i) 
  //  y2[i] = (std::fabs(b[i]) > delta) ? y[i] + a[i]/b[i]*(std::exp(b[i]*dt*0.5) - 1.0) : 
  //    y[i] + a[i]*dt*0.5;
  //
  //// Second step
  //
  //// Exact derivatives for linear terms
  //_ode->linear_derivatives(y2.data(), t, b.data());
  //_ode->eval(y2.data(), t, a.data());
  //
  //// Local variable to store comp i
  ////double yi;       
  //for (uint i = 0; i < num_states(); ++i) 
  //{        
  //  // Store original value of comp i
  //  //yi = y[i];
  //
  //  // Use y2[i] at diagonal
  //  //y[i] = y2[i];
  //
  //  // Evaluate right hand side after changing comp i    
  //  //a[i] = _ode->eval(i, y, t);
  //
  //  if (linear_terms[i] == 0) 
  //  {  
  //    // Numerical differentiation for non-linear terms
  //    y2[i] += delta; 
  //
  //    // Component i derivative
  //    b[i] = (_ode->eval(i, y2.data(), t) - a[i])/delta;
  //    y2[i] -= delta;                       // Restore state i
  //  }
  //
  //  // Restore state i
  //  //y[i] = yi;
  //}
  //
  //for (uint i = 0; i < num_states(); ++i) 
  //  y[i] += (std::fabs(b[i]) > delta) ? a[i]/b[i]*(std::exp(b[i]*dt) - 1.0) :
  //    a[i]*dt;
}
//-----------------------------------------------------------------------------
