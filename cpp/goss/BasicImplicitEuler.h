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

#ifndef BasicImplicitEuler_h_IS_INCLUDED
#define BasicImplicitEuler_h_IS_INCLUDED

#include <vector>
#include <memory>

#include "ImplicitODESolver.h"

namespace goss
{

  // Implicit Euler
  class BasicImplicitEuler : public ImplicitODESolver
  {
  public:

    // Default constructor
    BasicImplicitEuler();

    // Constructor
    BasicImplicitEuler(std::shared_ptr<ODE> ode);

    // Copy constructor
    BasicImplicitEuler(const BasicImplicitEuler& solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const
    { return std::make_shared<BasicImplicitEuler>(*this); }

    // Destructor
    ~BasicImplicitEuler ();

    // Attach ODE
    void attach(std::shared_ptr<ODE> ode);

    // Reset ODE
    void reset();

    // Solver specific compute jacobian method
    void compute_factorized_jacobian(double* y, double t, double dt);

    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);

  };

}
#endif
