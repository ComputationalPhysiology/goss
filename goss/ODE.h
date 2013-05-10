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

#ifndef ODE_H_IS_INCLUDED
#define ODE_H_IS_INCLUDED

#include <vector>
#include <boost/shared_ptr.hpp>

#include "types.h"
#include "DoubleVector.h"

namespace goss {

  // Base class for an ODE
  class ODE 
  {
  public:

    // Constructor
    ODE(uint num_states_);

    // Destructor
    virtual ~ODE() 
    {
      // Do nothing
    }

    // Return the size of the ODE
    inline uint num_states() const { return _num_states; }

    // Evaluate rhs of the ODE
    virtual void eval(const double* states, double t, double* values) = 0;

    // Evaluate component idx of the rhs of the ODE
    virtual double eval(uint idx, const double* states, double t);

    // Get default initial conditions
    virtual void get_ic(goss::DoubleVector* values) const = 0;

    // Return a copy of the ODE
    virtual boost::shared_ptr<ODE> copy() const = 0;

    // Compute numerical jacobian
    virtual void compute_jacobian(double t, double* states, double* jac);
    
    // In place LU Factorize matrix (jacobian)
    virtual void lu_factorize(double* mat);
    
    // Forward/Backward substitution of factoriesed matrix
    virtual void forward_backward_subst(const double* mat, const double* b, double* x);

  protected: 
    
    // ODE size
    const uint _num_states;

    
  private:

    // Temporaries used to compute jacobian
    std::vector<double> _f1, _f2;

  };
}

#endif
