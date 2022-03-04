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

#include "ESDIRK23a.h"
#include <math.h>
#include <stdio.h>

using namespace goss;

//-----------------------------------------------------------------------------
ESDIRK23a::ESDIRK23a()
    : AdaptiveImplicitSolver(), gamma(0.43586652150845899942), new_jacobian(true), first_step(true),
      a21(gamma), a22(gamma), a31((-4 * gamma * gamma + 6 * gamma - 1) / (4 * gamma)),
      a32((-2 * gamma + 1) / (4 * gamma)), a33(gamma), a41((6 * gamma - 1) / (12 * gamma)),
      a42(-1 / ((24 * gamma - 12) * gamma)),
      a43((-6 * gamma * gamma + 6 * gamma - 1) / (6 * gamma - 3)), a44(gamma), b1(a41), b2(a42),
      b3(a43), b4(a44), bh1(a31), bh2(a32), bh3(a33), c2(2.0 * gamma), c3(1.0), c4(1.0), s_fac(0.9),
      tol(1.0e-5), loc_error(1e9), prev_dt(0.1), eta(1.0), etamin(1.0e-10), kappa(1.0e-1),
      totalits(0), z1(0), z2(0), z3(0), z4(0), yn(0), yh(0), j_fac(0)
{
    _iord = 3;
    first_step = true;
}
//-----------------------------------------------------------------------------
ESDIRK23a::ESDIRK23a(std::shared_ptr<ODE> ode)
    : AdaptiveImplicitSolver(), gamma(0.43586652150845899942), new_jacobian(true), first_step(true),
      a21(gamma), a22(gamma), a31((-4 * gamma * gamma + 6 * gamma - 1) / (4 * gamma)),
      a32((-2 * gamma + 1) / (4 * gamma)), a33(gamma), a41((6 * gamma - 1) / (12 * gamma)),
      a42(-1 / ((24 * gamma - 12) * gamma)),
      a43((-6 * gamma * gamma + 6 * gamma - 1) / (6 * gamma - 3)), a44(gamma), b1(a41), b2(a42),
      b3(a43), b4(a44), bh1(a31), bh2(a32), bh3(a33), c2(2.0 * gamma), c3(1.0), c4(1.0), s_fac(0.9),
      tol(1.0e-5), loc_error(1e9), prev_dt(0.1), eta(1.0), etamin(1.0e-10), kappa(1.0e-1),
      totalits(0), z1(0), z2(0), z3(0), z4(0), yn(0), yh(0), j_fac(0)
{
    attach(ode);
    _iord = 3;
    first_step = true;
}
//-----------------------------------------------------------------------------
ESDIRK23a::ESDIRK23a(const ESDIRK23a &solver)
    : AdaptiveImplicitSolver(solver), gamma(0.43586652150845899942), new_jacobian(true),
      first_step(true), a21(gamma), a22(gamma),
      a31((-4 * gamma * gamma + 6 * gamma - 1) / (4 * gamma)), a32((-2 * gamma + 1) / (4 * gamma)),
      a33(gamma), a41((6 * gamma - 1) / (12 * gamma)), a42(-1 / ((24 * gamma - 12) * gamma)),
      a43((-6 * gamma * gamma + 6 * gamma - 1) / (6 * gamma - 3)), a44(gamma), b1(a41), b2(a42),
      b3(a43), b4(a44), bh1(a31), bh2(a32), bh3(a33), c2(2.0 * gamma), c3(1.0), c4(1.0), s_fac(0.9),
      tol(1.0e-5), loc_error(1e9), prev_dt(0.1), eta(1.0), etamin(1.0e-10), kappa(1.0e-1),
      totalits(0), z1(solver.num_states()), z2(solver.num_states()), z3(solver.num_states()),
      z4(solver.num_states()), yn(solver.num_states()), yh(solver.num_states()),
      j_fac(solver.num_states() * solver.num_states())
{
    _iord = 3;
}
//-----------------------------------------------------------------------------
ESDIRK23a::~ESDIRK23a()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void ESDIRK23a::attach(std::shared_ptr<ODE> ode)
{

    // Use base classes to actually attach ode
    ImplicitODESolver::attach(ode);

    z1.resize(num_states());
    z2.resize(num_states());
    z3.resize(num_states());
    z4.resize(num_states());
    yn.resize(num_states());
    yh.resize(num_states());
}
//-----------------------------------------------------------------------------
void ESDIRK23a::reset()
{
    // Reset counters
    nfevals = 0;
    ndtsa = 0;
    ndtsr = 0;

    // Reset bases
    AdaptiveImplicitSolver::reset();
}
//-----------------------------------------------------------------------------
void ESDIRK23a::compute_factorized_jacobian(const double &dt)
{
    /*This function assumes that the Jacobian of the rhs is already stored
    in jac. The factorized RK jacobian (I-dt*gamma*J) is stored in j_fac*/

    //std::cout << "Calling factorize jacobian" << std::endl;

    // FIXME: Hake, I think memcopy is better here...
    j_fac = _jac;

    mult(-dt * a22, j_fac.data());
    add_mass_matrix(j_fac.data());
    _ode->lu_factorize(j_fac.data());
    //recompute_jacobian = false;
}
//-----------------------------------------------------------------------------
void ESDIRK23a::advance_one_step(double *y, const double &t0, double &dt)
{

    bool accepted = do_step(y, t0, dt);
    while (!accepted && dt > min_dt) {
        if (!new_jacobian) {
            compute_ode_jacobian(y, t0); //, &jac[0]);
                                         //new_jacobian = true;
        }
        reduce_time_step(dt);
        compute_factorized_jacobian(dt);

        accepted = do_step(y, t0, dt);
    }
    new_jacobian = false;
}

bool ESDIRK23a::do_step(double *y, const double &t0, const double &delta_t)
{
    /*This function computes the 4 stage derivatives of the EESDIRK23a method.
    The first one is explicit while the rest are computed by calls to
    compute_stage_val.
    TODO:
    (i) compute_stage_val should be replaced by the more general
    newton_solve from ImplicitODESolver.
    (ii) The current formulation solves for stage derivatives. To reduce
    roundoff errors it should be reformulated to solve for (local) stage
    values instead */


    bool step_ok;
    uint i;

    _ode->eval(y, t0, z1.data());
    nfevals += 1;

    for (i = 0; i < num_states(); i++) {
        z2[i] = z1[i];
        _prev[i] = z1[i] * a21; //y[i]+z1[i]*a21*delta_t;
    }
    step_ok = compute_stage_val(z2, y, _prev, t0 + c2 * delta_t, delta_t);

    if (step_ok) {
        for (i = 0; i < num_states(); i++) {
            z3[i] = z2[i];
            _prev[i] = z1[i] * bh1 + z2[i] * bh2; //y[i]+delta_t*(z1[i]*bh1+z2[i]*bh2);
        }

        step_ok = compute_stage_val(z3, y, _prev, t0 + c3 * delta_t, delta_t);
    }

    if (step_ok) {
        for (i = 0; i < num_states(); i++) {
            z4[i] = z3[i];
            _prev[i] = z1[i] * b1 + z2[i] * b2
                       + z3[i] * b3; //y[i]+delta_t*(z1[i]*b1+z2[i]*b2+z3[i]*b3);
        }

        step_ok = compute_stage_val(z4, y, _prev, t0 + c4 * delta_t, delta_t);
    }

    if (step_ok) {
        //set new value of y = y0+delta_t*(b1*z1+...gamma*k4)
        for (i = 0; i < num_states(); i++) {
            //y[i] = y[i] + delta_t*(b1*z1[i]+b2*z2[i]+b3*z3[i]+gamma*z4[i]);
            yh[i] = delta_t
                    * ((b1 - bh1) * z1[i] + (b2 - bh2) * z2[i] + (b3 - gamma) * z3[i]
                       + gamma * z4[i]);
        }
        _ode->forward_backward_subst(j_fac.data(), yh.data(), yh.data());

        loc_error = norm(yh.data());
        if (loc_error > tol) {
            num_rejected++;
            //std::cout << "Loc error = " << loc_error << std::endl;
            step_ok = false;
        } else {
            for (i = 0; i < num_states(); i++)
                y[i] = y[i] + delta_t * (b1 * z1[i] + b2 * z2[i] + b3 * z3[i] + gamma * z4[i]);
        }
    }
    return step_ok;
}
//-----------------------------------------------------------------------------
bool ESDIRK23a::compute_stage_val(std::vector<double> &z, double *y, const std::vector<double> prev,
                                  const double &t_s, const double &dt)
{
    /*Local version of newton_solve, should be replaced by calls to the more
    generic function. */

    _newton_iterations = 0;
    uint i;
    double theta;
    bool conv_ok = true;

    std::vector<double> y_int(num_states()), z_int(num_states()), dz(num_states()),
            prev_dz(num_states());

    eta = pow(std::max(eta, etamin), 0.8);
    do {
        for (i = 0; i < num_states(); i++)
            y_int[i] = y[i] + dt * prev[i] + gamma * dt * z[i];


        _ode->eval(y_int.data(), t_s, z_int.data());


        for (i = 0; i < num_states(); i++)
            _b[i] = -z[i] + z_int[i];

        prev_dz = dz;

        // FIXME: This method does not do anything!
        compute_factorized_jacobian(y, t_s, dt);
        _ode->forward_backward_subst(j_fac.data(), _b.data(), dz.data());

        if (_newton_iterations >= 1) {
            theta = norm(dz.data()) / norm(prev_dz.data());
            if (!(theta <= 1)) {
                num_rejected++;
                conv_ok = false;
                break;
            }
            eta = theta / (1.0 - theta);
        }
        if (_newton_iterations > max_iterations) {
            //std::cout << "Not converged in "
            //	       << max_iterations <<" iterations, reducing time step" << std::endl;
            num_rejected++;
            conv_ok = false;
            break;
        }

        for (i = 0; i < num_states(); i++)
            z[i] += dz[i];
        _newton_iterations++;

    } while (eta * norm(dz.data()) > kappa * tol);
    totalits += _newton_iterations;

    return conv_ok;
}
//-----------------------------------------------------------------------------
void ESDIRK23a::reduce_time_step(double &delta_t)
{

    /*Reduces the time step based on the most recent value of the
    local error. Parameters are tuned for this particular ESDIRK23a method.*/

    double reduction_factor = 1.0;
    double e = loc_error;
    do {
        reduction_factor *= 0.6; //seems to be the best choice
        e *= 0.1;
    } while (e > tol);


    //std::cout << "Goss Time step reduced from " << delta_t;
    delta_t = std::max(reduction_factor * delta_t, min_dt);
    //delta_t = max(0.5*delta_t,min_dt);
    //std::cout << " to " << delta_t <<std::endl;
}
//-----------------------------------------------------------------------------
void ESDIRK23a::forward(double *y, const double t0, const double interval)
{

    double local_t = 0;

    //not used at this point...
    reached_tend = false;

    if (first_step) {
        _dt = interval;
        prev_dt = _dt;
        compute_ode_jacobian(y, t0); //,&jac[0]);
    } else {
        _dt = std::max(s_fac * prev_dt * pow((tol / loc_error), 0.25), min_dt);
        _dt = std::min(_dt, interval);
    }

    while (local_t < interval) {

        if (_dt * 1.1 >= interval - local_t)
            _dt = interval - local_t;

        if (_newton_iterations >= 2)
            compute_ode_jacobian(y, t0 + local_t); //, &jac[0]);

        if (fabs((_dt - prev_dt) / prev_dt) > 0.2 || new_jacobian)
            compute_factorized_jacobian(_dt);

        advance_one_step(y, t0 + local_t, _dt);

        prev_dt = _dt;
        local_t += _dt;
        _dt = std::max(s_fac * prev_dt * pow((tol / loc_error), 0.25), min_dt);
    }

    first_step = false;
    //std::cout << "@Rejects goss " << num_rejected << std::endl;
}
//-----------------------------------------------------------------------------
