// Copyright (C) 2012 Johan Hake
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

#ifndef LINEARIZED_ODE_H_IS_INCLUDED
#define LINEARIZED_ODE_H_IS_INCLUDED

#include "types.h"
#include "ODE.h"

namespace goss {

  // Class which provides an interface for ODE Solvers which need 
  // linearized terms (GRLX)
  class LinearizedODE : public virtual ODE
  {
  public:
    
    LinearizedODE(uint system_size) : ODE(system_size)
    { 
      // Do nothing
    } 

    virtual ~LinearizedODE() 
    {
      // Do nothing
    }

    // Populate indices with information about the linear terms
    virtual void linear_terms(uint* indices) const = 0;

    // Evaluate the linear derivatives
    virtual void linear_derivatives(const double* x, double t, double* y) const = 0;
    
  };
}

#endif
