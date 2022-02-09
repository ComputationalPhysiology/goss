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

#ifndef RL_H_IS_INCLUDED
#define RL_H_IS_INCLUDED

#include <vector>
#include <cmath>
#include <memory>

#include "types.h"
#include "ODESolver.h"
#include "ODE.h"

namespace goss {

  // First order accurate Generalized Rush-Larsen ODE Solver
  class RL1: public ODESolver
  {
  public:

    // Default Constructor
    RL1();

    // Constructor
    RL1(std::shared_ptr<ODE> ode);

    // Copy constructor
    RL1(const RL1& solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const { return std::make_shared<RL1>(*this); }

    // Destructor
    ~RL1();

    // Attach ODE to solver
    virtual void attach(std::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    virtual void forward(double* y, double t, double dt);

  protected:

    // One step of the RL algorithm
    inline void _one_step(double* y2, const double* y, const double* y0,
			  const double t, const double dt);

  };

  //-----------------------------------------------------------------------------
  inline void RL1::_one_step(double* y2, const double* y, const double* y0,
			     const double t, const double dt)
  {

    // Evaluate full right hand side
    _ode->linearized_eval(y, t, _f1().data(), _f2().data(), true);

    for (uint i = 0; i < num_states(); ++i)
      y2[i] = _ode->linear_term(i) ? y0[i] + _f2()[i]/_f1()[i]*(std::exp(_f1()[i]*dt) - 1.0) : y0[i] + _f2()[i]*dt;
  }

  //-----------------------------------------------------------------------------

}
#endif
