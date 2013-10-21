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
#include "log.h"
#include "constants.h"

using namespace goss;

//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver()
  : ODESolver(), jac(0), _f1(0), _yz(0), _b(0), _dz(0), _prev(0),
    _newton_tol(1.e-5), eta(1e-10), _kappa(0.1), jac_size(0), stages(0), newtonits(0), 
    maxits(10), rejects(0), jac_comp(0), num_tsteps(0), min_dt(0.0), 
    recompute_jacobian(true), _absolute_tol(1.e-10), _max_relative_residual(1.e-3)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver(double ldt)
  : ODESolver(ldt), jac(0), _f1(0), _yz(0), _b(0), _dz(0), _prev(0),
    _newton_tol(1.e-5), eta(1e-10), _kappa(0.1), jac_size(0), stages(0), newtonits(0), 
    maxits(10), rejects(0), jac_comp(0), num_tsteps(0), min_dt(0.0), 
    recompute_jacobian(true), _absolute_tol(1.e-10), _max_relative_residual(1.e-3)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver(const ImplicitODESolver& solver)
  : ODESolver(solver), jac(0), _f1(0), _yz(0), _b(0), _dz(0), _prev(0),
    _newton_tol(solver._newton_tol), eta(solver.eta), _kappa(solver._kappa), 
    jac_size(solver.jac_size), stages(solver.stages), newtonits(solver.newtonits), 
    maxits(solver.maxits), rejects(solver.rejects), jac_comp(solver.jac_comp), 
    num_tsteps(solver.num_tsteps), min_dt(solver.min_dt), 
    recompute_jacobian(solver.recompute_jacobian), _absolute_tol(solver._absolute_tol), 
    _max_relative_residual(solver._max_relative_residual)
{
  // Initialize memory
  _b.resize(num_states());
  _dz.resize(num_states());
  _prev.resize(num_states());

  _yz.resize(num_states());
  _f1.resize(num_states());

  // Init jacobian
  jac_size = num_states()*num_states();
  jac.resize(jac_size);

}
//-----------------------------------------------------------------------------
ImplicitODESolver::~ImplicitODESolver ()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::attach(boost::shared_ptr<ODE> ode)
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
  jac_size = num_states()*num_states();
  jac.resize(jac_size);

}
//-----------------------------------------------------------------------------
void ImplicitODESolver::reset()
{
  
  // Reset counts
  rejects = 0;
  jac_comp = 0;

  // We need to compute the Jacobian the first time
  recompute_jacobian = true;

  // Newton tolerance
  _newton_tol = 1.e-5;

  // Reset eta
  eta = 1.e-10;

}

//-----------------------------------------------------------------------------
void ImplicitODESolver::mult(double scale, double* mat)
{
  for (uint i = 0; i < num_states(); ++i)
    for (uint j = 0; j < num_states(); ++j)
      mat[i*num_states()+j] *= scale;
}
//-----------------------------------------------------------------------------
bool ImplicitODESolver::newton_solve(double* z, double* prev, double* y0, double t,
				     double dt, double alpha)
{
  uint i;
  bool step_ok = true;
  newtonits = 0;
  double relative_residual = 1.0, residual, prev_residual = 1.0;
  recompute_jacobian = false;

  do
  {

    // Compute solution 
    for (i = 0; i < num_states(); ++i)
      _yz[i] = y0[i] + z[i];

    // Evaluate ODE using computed solution
    _ode->eval(&_yz[0], t, &_f1[0]);
    
    // Build rhs for linear solve
    for (i = 0; i < num_states(); ++i)
      _b[i] = -z[i]*static_cast<double>(_ode->differential_states()[i]) + \
	dt*(prev[i] + alpha*_f1[i]);

    // Linear solve on factorized jacobian
    _ode->forward_backward_subst(&jac[0], &_b[0], &_dz[0]);
    residual = norm(&_dz[0]);

    // Check for residual convergence
    if (residual < _absolute_tol)
    {
      break;
    }

    // 2nd time around
    if (newtonits > 0) 
    {

      // How fast are we converging?
      relative_residual = residual/prev_residual;

      // If too slow we flag the jacobian to be recomputed
      recompute_jacobian = relative_residual >= _max_relative_residual;

      // If we diverge
      if (relative_residual > 1)
      {
	goss_debug1("Newton solver diverges with relative_residual: %f. "\
		    "Reducing time step.", relative_residual);
        rejects ++;
        step_ok = false;
	recompute_jacobian = true;
        break;
      }
      
      // We converge too slow
      if (residual > (_kappa*_newton_tol*(1 - relative_residual)/std::pow(relative_residual, maxits - newtonits)))
      {
	goss_debug3("Newton solver converges to slow with relative_residual: "\
		    "%.2e and residual: %.2e at iteration %d. "\
		    "Reducing time step.", relative_residual, residual, newtonits);
        rejects ++;
        step_ok = false;
        recompute_jacobian = true;
        break;
      }
      
      eta = relative_residual/(1.0 - relative_residual);

    }
    
    // newtonits == 0
    else
    {
      // On first iteration we need an approximation of eta. We take
      // the one from previous step and increase it slightly. This is
      // important for linear problems which only should recquire 1
      // iteration to converge.
      eta = eta > GOSS_EPS ? eta : GOSS_EPS;
      eta = std::pow(eta, 0.8);
    }

    // No convergence
    if (newtonits > maxits)
    {
      goss_debug1("Newton solver did not converged in %d iterations. Reducing " \
		  "time step.", maxits);
      recompute_jacobian = true;
      rejects ++;
      step_ok = false;
      //return step_ok;
      break;
    }
    
    // Update solution
    for (i = 0; i <num_states(); ++i)
      z[i] += _dz[i];

    prev_residual = residual;
    newtonits++;
    
    // eta*residul is the iteration error and an estimation of the
    // local discretization error.
  } while(eta*residual >= _kappa*_newton_tol);

  //goss_debug1("Newton converged in %d iterations.", newtonits);

  return step_ok;
}
//-----------------------------------------------------------------------------
double ImplicitODESolver::norm(double* vec)
{
  double l2_norm = 0;

  for (uint i = 0; i < num_states(); ++i)
    l2_norm += vec[i]*vec[i];

  l2_norm = std::sqrt(l2_norm);
  return l2_norm;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::add_mass_matrix(double* mat) const
{
  for (uint i = 0; i < num_states(); ++i)
    mat[i*num_states()+i] += static_cast<double>(_ode->differential_states()[i]);
}
//-----------------------------------------------------------------------------
