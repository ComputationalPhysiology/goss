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

#ifndef RL2_H_IS_INCLUDED
#define RL2_H_IS_INCLUDED

#include <vector>
#include <memory>

#include "types.h"
#include "ODESolver.h"
#include "RL1.h"

namespace goss {

  // Second order accurate Generalized Rush-Larsen ODE Solver
  class RL2 : public RL1
  {
  public:

    // Default Constructor
    RL2();

    // Constructor
    RL2(std::shared_ptr<ODE> ode);

    // Copy constructor
    RL2(const RL2& solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const { return std::make_shared<RL2>(*this); }

    // Destructor
    ~RL2();

    // Attach ODE to solver
    virtual void attach(std::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);

  private:

    // Pointers to intermediate values used while stepping
    std::vector<double> _y2;

  };

}
#endif
