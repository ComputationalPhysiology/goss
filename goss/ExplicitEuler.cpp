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

#include "ExplicitEuler.h"

using namespace goss;

//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler() : ODESolver(), _dFdt(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(double ldt) : ODESolver(ldt), _dFdt(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(boost::shared_ptr<ODE> ode, double ldt) : 
  ODESolver(ldt), _dFdt(0)
{ 
  attach(ode);
} 
//-----------------------------------------------------------------------------
ExplicitEuler::ExplicitEuler(const ExplicitEuler& solver) : 
  ODESolver(solver), _dFdt(new double[solver.num_states()])
{ 
  // Do nothing
} 
//-----------------------------------------------------------------------------
ExplicitEuler::~ExplicitEuler()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ExplicitEuler::attach(boost::shared_ptr<ODE> ode)
{
  // Attach ODE
  ODESolver::attach(ode);

  // Create memory for derivative evaluation
  _dFdt.reset(new double[num_states()]);

}
//-----------------------------------------------------------------------------
void ExplicitEuler::forward(double* y, double t, double interval) 
{

  assert(_ode);
  
  // Calculate number of steps and size of timestep based on _ldt
  const ulong nsteps = _ldt > 0 ? std::ceil(interval/_ldt - 1.0E-12) : 1;
  const double dt = interval/nsteps;

  // Local time
  double lt = t;
  for (ulong step = 0; step < nsteps; ++step)
  {
    // Evaluate rhs
    _ode->eval(y, lt, _dFdt.get());

    // Update states
    for (uint i = 0; i < num_states(); ++i)
      y[i] += dt*_dFdt[i];

    // Increase time
    lt += dt;

  }

}
