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

#ifndef ESDIRK23a_h_IS_INCLUDED
#define ESDIRK23a_h_IS_INCLUDED

#include <memory>
#include <vector>

#include "AdaptiveImplicitSolver.h"

namespace goss
{

// Explicit Singly Diagonally Implicit Runge-Kutta solver
class ESDIRK23a : public AdaptiveImplicitSolver
{

  public:
    // Default parameters
    int num_refinements_without_always_recomputing_jacobian = 2;
    double min_dt = 0.001;

    //Default constructor
    ESDIRK23a();

    // Constructor
    ESDIRK23a(std::shared_ptr<ODE> ode);

    // Copy constructor
    ESDIRK23a(const ESDIRK23a &solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const
    {
        return std::make_shared<ESDIRK23a>(*this);
    }

    // Attach ODE
    virtual void attach(std::shared_ptr<ODE> ode);

    // Jacobian of rhs
    virtual void compute_ode_jacobian(double *y, double t)
    {
        _ode->compute_jacobian(y, t, _jac.data());
        _jac_comp++;
        new_jacobian = true;
    }

    //virtual void compute_factorized_jacobian(double* y, double t, double dt){}
    virtual void compute_factorized_jacobian(const double &dt);

    //not really needed:
    virtual void compute_factorized_jacobian(double *, double, double){};

    // Reset ODE
    virtual void reset();

    // Step solver an interval of time forward
    void forward(double *y, double t, double interval);

    // Destructor
    ~ESDIRK23a();

    // Counters for the number of right hand side evaluations (nfevals) and
    // the number of accepted and rejected timesteps (ndtsa, ndtsr)
    long nfevals, ndtsa, ndtsr;

  private:
    void advance_one_step(double *y, const double &t0, double &dt);
    bool do_step(double *y, const double &t0, const double &delta_t);
    bool compute_stage_val(std::vector<double> &z, double *y, const std::vector<double> prev_stages,
                           const double &t_s, const double &dt);


    void reduce_time_step(double &delta_t);
    // Help variable
    double gamma;

    bool new_jacobian, first_step;

    // RK coefficients
    double a21, a22, a31, a32, a33, a41, a42, a43, a44;

    // RK weights
    double b1, b2, b3, b4, bh1, bh2, bh3;

    // RK coefficients
    double c2, c3, c4;

    //used in adaptive time step algorithm
    double s_fac, tol, loc_error, prev_dt; //, min_dt;
    double eta, etamin, kappa;
    uint totalits; //work, its, maxits;
    // State derivatives, allocated in attach(ode)
    std::vector<double> z1, z2, z3, z4, yn, yh;

    //vector to store the factorized matrix (I-dt*gamma*J)
    std::vector<double> j_fac;
};

} // namespace goss
#endif
