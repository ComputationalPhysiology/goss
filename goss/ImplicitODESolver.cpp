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

// FIXME: This class needs some serious reconsideration regarding memory use!!

#include "ImplicitODESolver.h"
#include <iostream>
#include <cmath>

using namespace goss;

//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver()
  : ODESolver(), jac(0), _f1(0), _f2(0), _yz(0), _b(0), _dz(0), _prev(0),
    _newton_tol(1.e-5), eta(1e-10), kappa(0.1), stages(0), newtonits(0), 
    maxits(10), rejects(0), jac_comp(0), min_dt(0.0), recompute_jacobian(true)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver(double ldt)
  : ODESolver(ldt), jac(0), _f1(0), _f2(0), _yz(0), _b(0), _dz(0), _prev(0),
    _newton_tol(1.e-5), eta(1e-10), kappa(0.1), stages(0), newtonits(0), 
    maxits(10), rejects(0), jac_comp(0), min_dt(0.0), recompute_jacobian(true)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver(const ImplicitODESolver& solver)
  : ODESolver(solver), jac(0), _f1(0), _f2(0), _yz(0), _b(0), _dz(0), _prev(0),
    _newton_tol(solver._newton_tol), eta(solver.eta), kappa(solver.kappa), 
    stages(solver.stages), newtonits(solver.newtonits), 
    maxits(solver.maxits), rejects(solver.rejects), jac_comp(solver.jac_comp), 
    min_dt(solver.min_dt), recompute_jacobian(solver.recompute_jacobian)
{
  // Initialize memory
  _b.resize(num_states());
  _dz.resize(num_states());
  _prev.resize(num_states());

  _yz.resize(num_states());
  _f1.resize(num_states());
  _f2.resize(num_states());

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
  _f2.resize(num_states());

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

}

//-----------------------------------------------------------------------------
void ImplicitODESolver::mult(double scale, double* mat)
{
  for (uint i = 0; i < num_states(); ++i)
    for (uint j = 0; j < num_states(); ++j)
      mat[i*num_states()+j] *= scale;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::add_identity(double* mat)
{
  for (uint i = 0; i < num_states(); ++i)
    mat[i*num_states()+i] += 1;
}
//-----------------------------------------------------------------------------
bool ImplicitODESolver::newton_solve(double* z, double* prev, double* y0, double t,
				     double dt, double alpha)
{
  uint i;
  bool step_ok = true, converged = false;
  newtonits = 0;
  double Ntheta, z_norm, prev_norm = 1.0;
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
      _b[i] = -z[i] + dt*(prev[i] + alpha*_f1[i]);

    // Linear solve on factorized jacobian
    _ode->forward_backward_subst(&jac[0], &_b[0], &_dz[0]);
    z_norm = norm(&_dz[0]);

    // 2nd time around
    if (newtonits > 0) 
    {

      // How fast are we converging?
      Ntheta = z_norm/prev_norm;

      // If not fast enough recompute jacobian
      if (Ntheta < 1e-3)
        recompute_jacobian = false;
      else
        recompute_jacobian = true;
    
      // If we diverge
      if (Ntheta > 1)
      {
#ifdef DEBUG
        std::cout << "Newton solver diverges with Ntheta = " << Ntheta << \
	  ", reduces time step." << std::endl;
#endif
        rejects ++;
        step_ok = false;
        //return step_ok;
        break;
      }
      
      // We converge too slow
      if (z_norm > (kappa*_newton_tol*(1 - Ntheta)/std::pow(Ntheta, maxits - newtonits)))
      {
#ifdef DEBUG
        std::cout << "Newton solver converges to slow with Ntheta = " << Ntheta << \
	  " at iteration " << newtonits << ", reduces time step." << std::endl;
#endif
        rejects ++;
        step_ok = false;
        recompute_jacobian = true;
        break;
      }
      
      eta = Ntheta/(1.0 - Ntheta);
    }
    
    // newtonits == 0
    else
    {
      eta = eta > 1e-15 ? eta : 1e-15;
      eta = std::pow(eta, 0.8);
    }

    // No convergence
    if (newtonits > maxits && !converged)
    {
#ifdef DEBUG
      std::cout << "Not converged in " << maxits << " iterations. Reduces "\
	"time step." << std::endl;
#endif
      rejects ++;
      step_ok = false;
      //return step_ok;
      break;
    }
    
    // Update solution
    for (i = 0; i <num_states(); ++i)
      z[i] += _dz[i];

    prev_norm = z_norm;
    newtonits++;

  } while (eta*z_norm <= kappa*_newton_tol); 

#ifdef DEBUG
  std::cout << "Newton converged in " << newtonits << " iterations." << std::endl;
#endif
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
