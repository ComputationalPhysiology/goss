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

#include <memory>
#include <vector>

#include "DoubleVector.h"
#include "types.h"

namespace goss
{

// Base class for an ODE
class ODE
{
  public:
      // Constructor
    ODE(uint num_states_);

    // Copy constructor
    ODE(const ODE &ode);

    // Destructor
    virtual ~ODE()
    {
        // Do nothing
    }

    // Return the size of the ODE
    inline uint num_states() const
    {
        return _num_states;
    }

    // Evaluate rhs of the ODE
    virtual void eval(const double *states, double time, double *values) = 0;

    // Evaluate component idx of the rhs of the ODE
    virtual double eval(uint idx, const double *states, double time);

    // Evaluate the linearized rhs
    virtual void linearized_eval(const double *states, double time, double *linearized, double *rhs,
                                 bool only_linear) const;

    // Get default initial conditions
    virtual void get_ic(goss::DoubleVector *values) const = 0;

    // Return a copy of the ODE
    virtual std::shared_ptr<ODE> copy() const = 0;

    // Compute numerical jacobian
    virtual void compute_jacobian(double *states, double time, double *jac);

    // In place LU Factorize matrix (jacobian)
    virtual void lu_factorize(double *mat) const;

    // Forward/Backward substitution of factorized matrix
    virtual void forward_backward_subst(const double *mat, const double *b, double *x) const;

    // Returns true if the ODE is a DAE (then only implicit solvers can be used)
    bool is_dae() const
    {
        return _is_dae;
    }

    // Return a view of the differential states of the ODE
    const std::vector<unsigned char> &differential_states() const
    {
        return _differential_states;
    }

    // Return wether the ith term is linear
    inline unsigned char linear_term(uint i) const
    {
        return _linear_terms[i];
    }

  protected:
    // ODE size
    const uint _num_states;

    // Flags for what states are differential
    std::vector<unsigned char> _differential_states;

    // Flags for what states the rhs is linearly dependent on itself
    std::vector<unsigned char> _linear_terms;

    // Flag to determine if an ODE is a DAE
    bool _is_dae;

  private:
    // Friends
    friend class ODESolver;

    // Temporaries used to compute jacobian and other temporaries
    std::vector<double> _f1, _f2;
};
} // namespace goss

#endif
