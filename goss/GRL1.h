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

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "types.h"
#include "ODESolver.h"
#include "ODE.h"

namespace goss {

  // First order accurate Generalized Rush-Larsen ODE Solver
  class GRL1: public ODESolver
  {
  public:
    
    // Default Constructor
    GRL1();

    // Constructor
    GRL1(boost::shared_ptr<ODE> ode);
    
    // Copy constructor
    GRL1(const GRL1& solver);
    
    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const { return boost::make_shared<GRL1>(*this); }

    // Destructor
    ~GRL1();

    // Attach ODE to solver
    virtual void attach(boost::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);
    
  private:

    // Pointers to intermediate values used while stepping
    std::vector<double> a;
    std::vector<double> b;

    std::vector<uint> linear_terms;
    const double delta;

  };
}
#endif
