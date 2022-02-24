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
#include <cstdio>
#include <cstring>

#include "GRL2.h"
#include "Timer.h"
#include "log.h"

using namespace goss;

//-----------------------------------------------------------------------------
GRL2::GRL2() : GRL1(), _y2(0)
{
}
//-----------------------------------------------------------------------------
GRL2::GRL2(std::shared_ptr<ODE> ode) : GRL1(), _y2(0)
{
    attach(ode);
}
//-----------------------------------------------------------------------------
GRL2::GRL2(const GRL2 &solver) : GRL1(solver), _y2(solver._y2)
{
}
//-----------------------------------------------------------------------------
GRL2::~GRL2()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void GRL2::attach(std::shared_ptr<ODE> ode)
{
    // Attach ode using base class
    ODESolver::attach(ode);

    if (ode->is_dae())
        goss_error("GRL2.cpp", "attaching ode",
                   "cannot integrate a DAE ode with an explicit solver.");

    // Initalize memory
    _y2.resize(num_states(), 0.0);
}
//-----------------------------------------------------------------------------
void GRL2::forward(double *y, double t, double dt)
{
    // Calculate number of steps and size of timestep based on _ldt
    const double ldt_0 = _ldt;
    const double delta = delta;
    const ulong nsteps = ldt_0 > 0 ? std::ceil(dt / ldt_0 - 1.0E-12) : 1;
    const double ldt = dt / nsteps;

    // Local time
    double lt = t;

    for (ulong step = 0; step < nsteps; ++step) {

        // First step
        _one_step(_y2.data(), y, y, lt, ldt * 0.5, delta);

        // Second step
        _one_step(y, _y2.data(), y, lt, ldt, delta);

        // Increase time
        lt += ldt;
    }
}
//-----------------------------------------------------------------------------
