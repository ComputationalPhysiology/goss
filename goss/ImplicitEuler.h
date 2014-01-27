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

#ifndef ImplicitEuler_h_IS_INCLUDED
#define ImplicitEuler_h_IS_INCLUDED

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "ImplicitODESolver.h"

namespace goss 
{

  // Implicit Euler
  class ImplicitEuler : public ImplicitODESolver
  {
  public:

    // Default constructor
    ImplicitEuler();

    // Constructor
    ImplicitEuler(boost::shared_ptr<ODE> ode, double ldt=-1.0);

    // Constructor
    ImplicitEuler(double ldt);

    // Copy constructor
    ImplicitEuler(const ImplicitEuler& solver);

    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const 
    { return boost::make_shared<ImplicitEuler>(*this); }

    // Destructor
    ~ImplicitEuler ();

    // Attach ODE
    void attach(boost::shared_ptr<ODE> ode);

    // Reset ODE
    void reset();

    // Solver specific compute jacobian method
    void compute_factorized_jacobian(double* y, double t, double dt);

    // Step solver an interval of time forward
    void forward(double* y, double t, double interval);

    std::vector<int> newton_iter1;
    std::vector<int> newton_accepted1;
    std::vector<double> dt_v;

  protected:

    std::vector<double> z1;
    bool justrefined;

  };

}
#endif
