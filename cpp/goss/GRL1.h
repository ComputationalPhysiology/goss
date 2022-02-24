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

#ifndef GRL1_H_IS_INCLUDED
#define GRL1_H_IS_INCLUDED

#include <cmath>
#include <memory>
#include <vector>

#include "ODE.h"
#include "ODESolver.h"
#include "types.h"

namespace goss
{

// First order accurate Generalized Rush-Larsen ODE Solver
class GRL1 : public ODESolver
{
  public:
    // Default parameters
    double delta = 1e-8;

    // Default Constructor
    GRL1();

    // Constructor
    GRL1(std::shared_ptr<ODE> ode);

    // Copy constructor
    GRL1(const GRL1 &solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const
    {
        return std::make_shared<GRL1>(*this);
    }

    // Destructor
    ~GRL1();

    // Attach ODE to solver
    virtual void attach(std::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    virtual void forward(double *y, double t, double dt);

  protected:
    // One step of the GRL algorithm
    inline void _one_step(double *y2, const double *y, const double *y0, const double t,
                          const double dt, const double delta);
};

//-----------------------------------------------------------------------------
inline void GRL1::_one_step(double *y2, const double *y, const double *y0, const double t,
                            const double dt, const double delta)
{
    assert(_ode);

    // Evaluate full right hand side
    _ode->linearized_eval(y, t, _f1().data(), _f2().data(), false);

    for (uint i = 0; i < num_states(); ++i)
        y2[i] = (std::fabs(_f1()[i]) > delta)
                        ? y0[i] + _f2()[i] / _f1()[i] * (std::exp(_f1()[i] * dt) - 1.0)
                        : y0[i] + _f2()[i] * dt;
}

//-----------------------------------------------------------------------------

} // namespace goss
#endif
