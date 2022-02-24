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

#include <cmath>

#include "ImplicitODESolver.h"
#include "constants.h"
#include "log.h"

using namespace goss;

//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver()
    : ODESolver(), _jac(0), _f1(0), _yz(0), _b(0), _dz(0), _prev(0), _eta(1.), _stages(0),
      _rejects(0), _jac_comp(0), _recompute_jacobian(true), _newton_iterations(0)
{
}
//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver(const ImplicitODESolver &solver)
    : ODESolver(solver), _jac(0), _f1(0), _yz(0), _b(0), _dz(0), _prev(0), _eta(solver._eta),
      _stages(solver._stages), _rejects(solver._rejects), _jac_comp(solver._jac_comp),
      _recompute_jacobian(solver._recompute_jacobian), _newton_iterations(solver._newton_iterations)
{

    // Initialize memory
    _b.resize(num_states());
    _dz.resize(num_states());
    _prev.resize(num_states());

    _yz.resize(num_states());
    _f1.resize(num_states());

    // Init jacobian
    _jac.resize(num_states() * num_states());
}
//-----------------------------------------------------------------------------
ImplicitODESolver::~ImplicitODESolver()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::attach(std::shared_ptr<ODE> ode)
{

    // Attach ode using base attach
    // NOTE: This calls reset in the most derived class, which then propagates
    // NOTE: the call downwards to base classes.
    ODESolver::attach(ode);

    // Initialize memory
    _b.resize(num_states());
    _dz.resize(num_states());
    _prev.resize(num_states());

    _yz.resize(num_states());
    _f1.resize(num_states());

    // Init jacobian
    _jac.resize(num_states() * num_states());
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::reset()
{

    // Reset counts
    _rejects = 0;
    _jac_comp = 0;

    // We need to compute the Jacobian the first time
    _recompute_jacobian = true;

    // Reset eta
    _eta = eta_0;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::compute_factorized_jacobian(double *y, double t, double dt, double alpha)
{

    // Let ODE compute the jacobian
    _ode->compute_jacobian(y, t, _jac.data());

    // Build scaled discretization of jacobian
    mult(-dt * alpha, _jac.data());
    add_mass_matrix(_jac.data());

    // Factorize the jacobian
    _ode->lu_factorize(_jac.data());
    _jac_comp += 1;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::mult(double scale, double *mat)
{
    for (uint i = 0; i < num_states(); ++i)
        for (uint j = 0; j < num_states(); ++j)
            mat[i * num_states() + j] *= scale;
}
//-----------------------------------------------------------------------------
bool ImplicitODESolver::newton_solve(double *z, double *prev, double *y0, double t, double dt,
                                     double alpha, bool always_recompute_jacobian)
{
    uint i;
    bool step_ok = true;
    _newton_iterations = 0;
    double relative_previous_residual = 1.0;
    double relative_residual = 1.0;
    double previous_residual = 1.0;
    double initial_residual = 1.0;
    double residual;

    do {

        // Compute local newton solution
        for (i = 0; i < num_states(); ++i)
            _yz[i] = y0[i] + z[i];

        // Evaluate ODE using local solution
        _ode->eval(_yz.data(), t, _f1.data());

        // Build rhs for linear solve
        // z = y-y0
        // prev is a linear combination of previous stage solutions
        for (i = 0; i < num_states(); ++i)
            _b[i] = -z[i] * _ode->differential_states()[i] + dt * (prev[i] + alpha * _f1[i]);

        // Calculate the residual
        residual = norm(_b.data());

        // Check for relative residual convergence
        if (relative_residual < relative_tolerance)
            break;

        // Recompute jacobian if nessecary
        if (_recompute_jacobian || always_recompute_jacobian) {
            compute_factorized_jacobian(_yz.data(), t, dt, alpha);
            _recompute_jacobian = false;
        }

        // Linear solve on factorized jacobian
        _ode->forward_backward_subst(_jac.data(), _b.data(), _dz.data());

        // _Newton_Iterations == 0
        if (_newton_iterations == 0) {

            initial_residual = residual;

            // On first iteration we need an approximation of eta. We take
            // the one from previous step and increase it slightly. This is
            // important for linear problems which only should recquire 1
            // iteration to converge.
            _eta = _eta > GOSS_EPS ? _eta : GOSS_EPS;
            _eta = std::pow(_eta, 0.8);
        }

        // 2nd time around
        else {

            // How fast are we converging?
            relative_previous_residual = residual / previous_residual;

            // If too slow we flag the jacobian to be recomputed
            _recompute_jacobian = relative_previous_residual >= max_relative_previous_residual;

            // If we diverge
            if (relative_previous_residual >= 1) {
                log(DBG,
                    "Diverges       | t : %g, it : %2d, relative_previous_residual: %f, "
                    "relativ_residual: %g. Reducing time step and recompute jacobian.",
                    t, _newton_iterations, relative_previous_residual, relative_residual);
                step_ok = false;
                _rejects++;
                _recompute_jacobian = true;
                break;
            }

            const double scaled_relative_previous_residual = std::max(
                    std::pow(relative_previous_residual, max_iterations - _newton_iterations),
                    GOSS_EPS);
            // We converge too slow
            if (residual > (kappa * relative_tolerance * (1 - relative_previous_residual)
                            / scaled_relative_previous_residual)) {
                log(DBG,
                    "To slow        | t : %g, it: %2d, relative_previous_residual: "
                    "%f, relative_residual: %g. Recomputing Jacobian.",
                    t, _newton_iterations, relative_previous_residual, relative_residual);
                _recompute_jacobian = true;
            }

            _eta = relative_previous_residual / (1.0 - relative_previous_residual);
        }

        // No convergence
        if (_newton_iterations > max_iterations) {
            log(DBG,
                "Max iterations | t : %g, it: %2d, relative_previous_residual: "
                "%f, relative_residual: %g. Recomputing Jacobian.",
                t, _newton_iterations, relative_previous_residual, relative_residual);
            _recompute_jacobian = true;
            _rejects++;
            step_ok = false;
            break;
        }

        // Update local increment solution
        for (i = 0; i < num_states(); ++i)
            z[i] += _dz[i];

        // Update residuals
        relative_residual = residual / initial_residual;
        previous_residual = residual;
        _newton_iterations++;

        log(5,
            "Monitor        | t : %g, it : %2d, relative_previous_residual: %f, "
            "relativ_residual: %g.",
            t, _newton_iterations, relative_previous_residual, relative_residual);
        // eta*residul is the iteration error and an estimation of the
        // local discretization error.
    } while (_eta * relative_residual >= kappa * relative_tolerance);

    //goss_debug1("Newton converged in %d iterations.", _newton_iterations);

    return step_ok;
}
//-----------------------------------------------------------------------------
double ImplicitODESolver::norm(double *vec)
{
    double l2_norm = 0;

    for (uint i = 0; i < num_states(); ++i)
        l2_norm += vec[i] * vec[i];

    l2_norm = std::sqrt(l2_norm);
    return l2_norm;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::add_mass_matrix(double *mat, double weight) const
{
    for (uint i = 0; i < num_states(); ++i)
        mat[i * num_states() + i] += weight * _ode->differential_states()[i];
}
//-----------------------------------------------------------------------------
