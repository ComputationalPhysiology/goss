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

#ifndef ESDIRK4O32_h_IS_INCLUDED
#define ESDIRK4O32_h_IS_INCLUDED

#include <vector>
#include <memory>

#include "AdaptiveImplicitSolver.h"

namespace goss
{

  // Explicit Singly Diagonally Implicit Runge-Kutta solver
  class ESDIRK4O32: public AdaptiveImplicitSolver
  {

    public:

    // Default constructor
    ESDIRK4O32();

    // Constructor
    ESDIRK4O32(std::shared_ptr<ODE> ode);

    // Copy constructor
    ESDIRK4O32(const ESDIRK4O32& solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const
    { return std::make_shared<ESDIRK4O32>(*this); }

    // Attach ODE
    virtual void attach(std::shared_ptr<ODE> ode);

    // Reset ODE
    virtual void reset();

    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);

    // Destructor
    ~ESDIRK4O32 ();

    // Counters for the number of right hand side evaluations (nfevals) and
    // the number of accepted and rejected timesteps (ndtsa, ndtsr)
    long nfevals, ndtsa, ndtsr;

  private:

    // Help variable
    double gamma;

    // RK coefficients
    double a21, a22, a31, a32, a33, a41, a42, a43, a44;

    // RK weights
    double b1, b2, b3, b4, bh1, bh2, bh3;

    // RK coefficients
    double c2, c3, c4;

    // State derivatives, allocated in attach(ode)
    std::vector<double> z1, z2, z3, z4, yn, yh;

  };

}
#endif
