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

#ifndef RKF32_h_IS_INCLUDED
#define RKF32_h_IS_INCLUDED

#include <vector>
#include <memory>

#include "AdaptiveExplicitSolver.h"

namespace goss
{

  // Adaptive and explicit RungeKutta Solver
  class RKF32 : public AdaptiveExplicitSolver
  {
  public:

    // Constructor
    RKF32();

    // Constructor
    RKF32 (std::shared_ptr<ODE> ode);

    // Copy constructor
    RKF32(const RKF32& solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const { return std::make_shared<RKF32>(*this); }

    // Constructor
    virtual ~RKF32();

    // Attach ODE
    virtual void attach(std::shared_ptr<ODE> ode);

    // Reset solver
    virtual void reset();

    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);

    // FIXME: Where is this used!?
    // Store timestep and accepted timestep
    void log_data(double dt, bool accepted);

    // Return a vector of collected timesteps
    void dt_vector(DoubleVector *res);

    // Return a record of accepted
    void accepted_vector(DoubleVector *res);

    // Counters for the number of right hand side evaluations (nfevals) and
    // the number of accepted and rejected timesteps (ndtsa, ndtsr)
    long nfevals, ndtsa, ndtsr;

  private:

    // RK coefficients
    const double a21, a32;

    // RK weights
    const double b1, b2, b3, bh1, bh2, bh3, bh4;

    // Error weights
    const double d1, d2, d3, d4;

    // RK nodes
    const double c2, c3;

    // System size in bytes
    ulong nbytes;

    // State derivatives, allocated in attach(ode)
    std::vector<double> ki, k1, k2, k3, k4, yn, e;

    // Parameter for scalar or vector tolerance computing
    bool first;
  };

}
#endif
