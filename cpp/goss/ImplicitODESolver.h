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

#include "ODESolver.h"
#include <memory>
#include <vector>

namespace goss
{

// Base class of all Implicit Solvers
class ImplicitODESolver : public ODESolver
{

  public:
    // Default parameters
    double eta_0 = 1.0;
    double kappa = 0.1;
    double relative_tolerance = 1.e-12;
    int max_iterations = 30;
    double max_relative_previous_residual = 0.01;
    bool always_recompute_jacobian = false;

    // Default constructor
    ImplicitODESolver();

    // Copy constructor
    ImplicitODESolver(const ImplicitODESolver &solver);

    // Destructor
    virtual ~ImplicitODESolver();

    // Reset solver
    virtual void reset();

    // Attach ode
    virtual void attach(std::shared_ptr<ODE> ode);

    // Step solver an interval of time forward
    virtual void forward(double *y, double t, double interval) = 0;

    // Solver specific compute jacobian method
    void compute_factorized_jacobian(double *y, double t, double dt, double alpha);

    // Return the number of recomputation of the jacobian
    int num_jac_comp()
    {
        return _jac_comp;
    }

  protected:
    // Scale a matrix
    void mult(double fact, double *mat);

    // Add the mass matrix based on what states are differential
    void add_mass_matrix(double *mat, double weight = 1.0) const;

    // This function is designed for SDIRK and Backward Euler:
    virtual bool newton_solve(double *k, double *prev, double *y0, double t, double dt,
                              double alpha, bool always_recompute_jacobian);

    // Compute the norm of a vector
    virtual double norm(double *vec);

    // Variables used in the jacobian evaluation
    std::vector<double> _jac;
    std::vector<double> _f1, _yz;

    // Right hand side and solution of the linear system
    std::vector<double> _b, _dz;

    // Previous stages, used by DIRK methods
    std::vector<double> _prev;

    // Variable used in the estimation of the error of the newton
    // iteration for the first iteration (Important for linear problems!)
    double _eta;

    uint _stages;
    int _rejects, _jac_comp;
    bool _recompute_jacobian;
    int _newton_iterations;
};

} // namespace goss
#endif
