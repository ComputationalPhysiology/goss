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

#include "RK4.h"
#include "log.h"

using namespace goss;

//-----------------------------------------------------------------------------
RK4::RK4() : ODESolver(), k1(0), k2(0), k3(0), k4(0), tmp(0)
{
}
//-----------------------------------------------------------------------------
RK4::RK4(std::shared_ptr<ODE> ode) : ODESolver(), k1(0), k2(0), k3(0), k4(0), tmp(0)
{
    attach(ode);
}
//-----------------------------------------------------------------------------
RK4::RK4(const RK4 &solver)
    : ODESolver(solver), k1(solver.num_states()), k2(solver.num_states()), k3(solver.num_states()),
      k4(solver.num_states()), tmp(solver.num_states())
{
    // Do nothing
}
//-----------------------------------------------------------------------------
RK4::~RK4()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void RK4::attach(std::shared_ptr<ODE> ode)
{

    // Attach ode using base class
    ODESolver::attach(ode);

    if (ode->is_dae())
        goss_error("RK4.cpp", "attaching ode",
                   "cannot integrate a DAE ode with an explicit solver.");

    k1.resize(num_states());
    k2.resize(num_states());
    k3.resize(num_states());
    k4.resize(num_states());
    tmp.resize(num_states());
}
//-----------------------------------------------------------------------------
void RK4::forward(double *y, double t, double dt)
{

    assert(_ode);

    // Calculate number of steps and size of timestep based on ldt
    const double ldt_0 = _ldt;
    const ulong nsteps = ldt_0 > 0 ? std::ceil(dt / ldt_0 - 1.0E-12) : 1;
    const double ldt = dt / nsteps;

    // Local time
    double lt = t;
    for (ulong j = 0; j < nsteps; ++j) {
        // Evaluate rhs and calculate intermediate derivatives
        _ode->eval(y, lt, &k1[0]);

        // Explicit Euler step
        axpy(&tmp[0], y, 0.5 * ldt, &k1[0]);

        _ode->eval(&tmp[0], lt + 0.5 * ldt, &k2[0]);
        axpy(&tmp[0], y, 0.5 * ldt, &k2[0]);

        _ode->eval(&tmp[0], lt + 0.5 * ldt, &k3[0]);
        axpy(&tmp[0], y, ldt, &k3[0]);

        _ode->eval(&tmp[0], lt + ldt, &k4[0]);

        for (uint i = 0; i < num_states(); ++i)
            y[i] += ldt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;

        // Update local time
        lt += ldt;
    }
}
//-----------------------------------------------------------------------------
