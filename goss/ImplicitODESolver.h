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

#ifndef ImplicitODESolver_h_IS_INCLUDED
#define ImplicitODESolver_h_IS_INCLUDED

#include <vector>
#include <boost/shared_ptr.hpp>
#include "ODESolver.h"

namespace goss 
{

  // Base class of all Implicit Solvers
  class ImplicitODESolver : public ODESolver 
  {
  
  public:
  
    // Default constructor
    ImplicitODESolver();
  
    // Constructor
    ImplicitODESolver(double ldt);

    // Copy constructor
    ImplicitODESolver(const ImplicitODESolver& solver);

    // Destructor
    virtual ~ImplicitODESolver ();

    // Reset solver
    virtual void reset();

    // Attach ode
    virtual void attach(boost::shared_ptr<ODE> ode);

    // Step solver an interval of time forward
    virtual void forward(double* y, double t, double interval) = 0;

    // Set newton tolerance
    void set_newton_tol(double newton_tol){ _newton_tol = newton_tol; }

    // Return the number of recomputation of the jacobian
    int num_jac_comp(){ return jac_comp; }

  protected: 

    // Scale a matrix
    void mult(double fact, double* mat);

    // Add identity to matrix
    void add_identity(double* mat);

    // This function is designed for SDIRK and Backward Euler:
    virtual bool newton_solve(double* k, double* prev, double* y0, double t, 
			      double dt, double alpha);

    // Compute the norm of a vector
    virtual double norm(double* vec);

    // Variables used in the jacobian evaluation
    std::vector<double> jac;
    std::vector<double> _f1, _f2, _yz;

    // Right hand side and solution of the linear system
    std::vector<double> _b, _dz;

    // Previous stages, used by DIRK methods
    std::vector<double> _prev;

    // Newton tolerance
    double _newton_tol;

    // Variable used in the estimation of the error of the newton 
    // iteration for the first iteration (Important for linear problems!)
    double eta;

    // The safety factor for the stopping criterion of the newton iteration
    double kappa;
    uint jac_size;

    uint stages;
    int newtonits, maxits, rejects, jac_comp;
    ulong num_tsteps;

    double min_dt;
    bool recompute_jacobian;

  };

}
#endif
