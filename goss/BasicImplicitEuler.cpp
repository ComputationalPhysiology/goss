#include <cassert>
#include <cmath>
#include <cstdio>

#include "BasicImplicitEuler.h"
#include "log.h"
#include "constants.h"
#include "fenv.h"

using namespace goss;

//-----------------------------------------------------------------------------
BasicImplicitEuler::BasicImplicitEuler() : ImplicitODESolver(), newton_iter1(0), 
				 newton_accepted1(0), dt_v(0), z1(0)
{ 
  // Do nothing
} 
//-----------------------------------------------------------------------------
BasicImplicitEuler::BasicImplicitEuler(boost::shared_ptr<ODE> ode, double ldt) : 
  ImplicitODESolver(ldt), newton_iter1(0), newton_accepted1(0), dt_v(0), 
  z1(0)
{ 
  attach(ode);
} 
//-----------------------------------------------------------------------------
BasicImplicitEuler::BasicImplicitEuler(double ldt) : ImplicitODESolver(ldt), newton_iter1(0), 
					   newton_accepted1(0), dt_v(0), z1(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BasicImplicitEuler::BasicImplicitEuler(const BasicImplicitEuler& solver) : 
  ImplicitODESolver(solver), newton_iter1(solver.newton_iter1), 
  newton_accepted1(solver.newton_accepted1), dt_v(solver.dt_v), 
  z1(solver.num_states())
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BasicImplicitEuler::~BasicImplicitEuler()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::attach(boost::shared_ptr<ODE> ode)
{

  // Attach ode using bases
  ImplicitODESolver::attach(ode);

  // Init memory
  z1.resize(num_states());

}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::reset()
{
  
  num_tsteps = 0;

  stages = 1;
  newton_iter1.clear();
  newton_accepted1.clear();
  dt_v.clear();
  
  ImplicitODESolver::reset();
}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::compute_factorized_jacobian(double* y, double t, double dt)
{
  
  // Let ODE compute the jacobian
  _ode->compute_jacobian(y, t, jac.data());

  // Build Euler discretization of jacobian
  mult(-dt, jac.data());
  add_mass_matrix(jac.data());

  // Factorize the jacobian
  _ode->lu_factorize(jac.data());
  jac_comp += 1;

}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::forward(double* y, double t, double interval)
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  assert(_ode);

  double residual, residual0=1.0, relative_residual; //incr_residual

  // Local time
  double lt = t+interval;
  double dt = interval;

  // NOT USED.
  bool recompute_jacobian_ = true;

  // Calculate number of steps and size of timestep based on _ldt

  /* USE THIS WHEN WE HAVE ONE NEWTON ITERATION INPLACE
  const ulong nsteps = _ldt > 0 ? std::ceil(interval/_ldt - 1.0E-12) : 1;
  const double dt = interval/nsteps;

  for (ulong step = 0; step < nsteps; ++step)
  {
    // Evaluate rhs
    _ode->eval(y, lt, &_dFdt[0]);

    // Update states
    for (uint i = 0; i < num_states(); ++i)
      y[i] += dt*_dFdt[i];

    // Increase time
    lt += dt;
  }

  for (i = 0; i < num_states(); ++i)
    _prev[i] = 0.0;
  */

  bool newton_converged = false;
  newtonits = 0;

  // Copy previous solution
  for (uint i = 0; i < num_states(); ++i)
    _prev[i] = y[i];
  
  // Evaluate ODE using computed solution
  _ode->eval(_prev.data(), lt, _f1.data());
    
  // Build rhs for linear solve
  // b = dt*eval(y) - M*(y-y0) 
  // With y=y0 at start of iteration
  // b = eval(y0)
  for (uint i = 0; i < num_states(); ++i)
    _b[i] = -dt*_f1[i];
  
  // Start iterations
  while (!newton_converged && newtonits < maxits)
  {
    
    // Compute Jacobian
    if (recompute_jacobian_)
      compute_factorized_jacobian(y, lt, dt);

    // Linear solve on factorized jacobian
    _ode->forward_backward_subst(jac.data(), _b.data(), _dz.data());
	  
    // Compute initial residual
    if (newtonits == 0)
    {
      residual = residual0 = norm(_b.data());
      //incr_residual = norm(_dz.data());
    }
    else
    {
      // Compute resdiual
      residual = norm(_b.data());
      //incr_residual = norm(_dz.data());
    }
	  
    for (uint i = 0; i < num_states(); ++i)
      y[i] -= _dz[i];
        
    // Update number of iterations
    ++newtonits;
        
    // Relative residual
    relative_residual = residual / residual0;

    // Output iteration number and residual
    log(DBG, "BasicImplicitEuler newton iteration %d: r (abs) = %.3e "	\
	"(tol = %.3e) r (rel) = %.3e (tol = %.3e)", newtonits,	\
	residual, _absolute_tol, relative_residual, _newton_tol);
    
    // Check convergence criteria
    if (relative_residual < _newton_tol || residual < _absolute_tol)
    {
      newton_converged = true;
    }
    else
    {
      // Evaluate ODE using computed solution
      _ode->eval(y, lt, _f1.data());
    
      // Build rhs for linear solve
      // b = eval(y) - (y-y0)/dt 
      for (uint i = 0; i < num_states(); ++i)
	_b[i] = (y[i]-_prev[i])*_ode->differential_states()[i] - dt*_f1[i];
    }
   
  }

  if (!newton_converged)
  {
    error("Newton solver did not converge. Maximal newton iterations exceded.");
  }
 
    
#ifdef DEBUG
  // Lower level than DEBUG!
  log(5, "BasicImplicitEuler done with comp_jac = %d at t=%1.2e\n", jac_comp, t);
#endif
}
//-----------------------------------------------------------------------------


