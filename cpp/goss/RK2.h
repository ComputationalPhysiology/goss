// Copyright (C) 2008-2012 Ola Skavhaug
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

#ifndef RK2_H_IS_INCLUDED
#define RK2_H_IS_INCLUDED

#include <vector>
#include <memory>

#include "ODESolver.h"
#include "types.h"

namespace goss
{

  // Explicit Runge Kutta solver of 2nd order
  class RK2 : public ODESolver
  {

  public:

    // Default constructor
    RK2();

    // Constructor
    RK2(std::shared_ptr<ODE> ode);

    // Copy constructor
    RK2(const RK2& solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const { return std::make_shared<RK2>(*this); }

    // Destructor
    ~RK2();

    // Attach ODE to solver
    void attach(std::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    virtual void forward(double* y, double t, double dt);

  protected:

    // State derivative, allocated in attach(ode)
    std::vector<double> k1, tmp;

    // Perform a weighted addition of y and z
    inline void axpy(double* x, const double* y, double a, const double* z);

  };
}
//-----------------------------------------------------------------------------
inline void goss::RK2::axpy(double* x, const double* y, double a, const double* z)
{
  for (uint i = 0; i < num_states(); ++i)
    x[i] = y[i] + a*z[i];
}
//-----------------------------------------------------------------------------
#endif
