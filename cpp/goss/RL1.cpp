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
#include <cstdlib>

#include "RL1.h"
#include "Timer.h"
#include "log.h"

using namespace goss;

//-----------------------------------------------------------------------------
RL1::RL1() : ODESolver()
{
}
//-----------------------------------------------------------------------------
RL1::RL1(std::shared_ptr<ODE> ode) : ODESolver()
{
    attach(ode);
}
//-----------------------------------------------------------------------------
RL1::RL1(const RL1 &solver) : ODESolver(solver)
{
}
//-----------------------------------------------------------------------------
RL1::~RL1()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void RL1::attach(std::shared_ptr<ODE> ode)
{
    // Attach ode using base class
    ODESolver::attach(ode);

    if (_ode->is_dae())
        goss_error("RL1.cpp", "attach ode", "cannot integrate a DAE ode with Rush Larsen method.");
}
//-----------------------------------------------------------------------------
void RL1::forward(double *y, double t, double dt)
{


    // Calculate number of steps and size of timestep based on _ldt
    const double ldt_0 = _ldt;
    const ulong nsteps = ldt_0 > 0 ? std::ceil(dt / ldt_0 - 1.0E-12) : 1;
    const double ldt = dt / nsteps;

    // Local time
    double lt = t;

    for (ulong step = 0; step < nsteps; ++step) {

        // One step
        _one_step(y, y, y, lt, ldt);

        // Increase time
        lt += ldt;
    }
}
//-----------------------------------------------------------------------------
