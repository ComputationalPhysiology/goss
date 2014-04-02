#include <stdio.h>
#include <cmath>

#include "ThetaSolver.h"
#include "log.h"

using namespace goss;

//-----------------------------------------------------------------------------
ThetaSolver::ThetaSolver() : ImplicitODESolver(), theta(0.5), 
			     newton_iter1(0), newton_accepted1(0), 
			     dt_v(0), z1(0), ft1(0), justrefined(false)
{}
//-----------------------------------------------------------------------------
ThetaSolver:: ThetaSolver (boost::shared_ptr<ODE> ode, double ldt) 
  : ImplicitODESolver(ldt), theta(0.5), newton_iter1(0), newton_accepted1(0), 
    dt_v(0), z1(0), ft1(0), justrefined(false)
{ 
  attach(ode);
} 
//-----------------------------------------------------------------------------
ThetaSolver::ThetaSolver(double ldt) : ImplicitODESolver(ldt), theta(0.5), 
				       newton_iter1(0), newton_accepted1(0), 
				       dt_v(0), z1(0), ft1(0), justrefined(false)
{}
//-----------------------------------------------------------------------------
void ThetaSolver::attach(boost::shared_ptr<ODE> ode)
{
  // Attach ode using bases
  ImplicitODESolver::attach(ode);

  // Init memory
  z1.resize(num_states());
  ft1.resize(num_states());

}
//-----------------------------------------------------------------------------
void ThetaSolver::reset()
{
  
  theta = 0.5;
  justrefined = false;
  num_tsteps = 0;

  stages = 1;
  newton_iter1.clear();
  newton_accepted1.clear();
  dt_v.clear();
  
  ImplicitODESolver::reset();
}
//-----------------------------------------------------------------------------
void ThetaSolver::compute_factorized_jacobian(double* y, double t, double dt)
{
  
  // Let ODE compute the jacobian
  _ode->compute_jacobian(y, t, &jac[0]);

  // Build Theta discretization of jacobian
  mult(-dt*(1-theta), &jac[0]);
  add_mass_matrix(&jac[0]);

  // Factorize the jacobian
  _ode->lu_factorize(&jac[0]);
  jac_comp += 1;

}
//-----------------------------------------------------------------------------
void ThetaSolver::forward(double* y, double t, double dt)
{

  assert(_ode);

  //std::cout << "theta: " << theta << " t: " << t << " dt: " << interval << " dt: " << _dt << " ldt: " << _ldt << " V: " << y[0] << std::endl;
  
  uint i;
  double t_end = t + dt;
  double ldt = _ldt > 0 ? _ldt : dt;
  
  for (i = 0; i < _ode->num_states(); ++i)
    _prev[i] = 0.0;

  bool step_ok, done = false;

  // A way to check if we are at t_end.
  const double eps = 1e-14;

  while (!done)
  {
    num_tsteps += 1;
  
    // Recompute the jacobian if nessesary
    if (recompute_jacobian)
    {
      compute_factorized_jacobian(y, t, ldt);
    }

    // Use 0.0 z1:
    for (i = 0; i < _ode->num_states(); ++i)
      z1[i] = 0.0;

    // Explicit eval
    _ode->eval(y, t, &ft1[0]);
    
    for (i = 0; i < _ode->num_states(); ++i)
      _prev[i] = theta*ft1[i];
    
    // Solve for increment
    step_ok = newton_solve(&z1[0], &_prev[0], y, t + ldt, ldt, theta);    
#ifdef DEBUG
    newton_iter1.push_back(newtonits);
    dt_v.push_back(ldt);
#endif

    // Newton step OK
    if (step_ok)
    {
    
      t+=ldt;
      
      if (std::fabs(t-t_end) < eps)
      {
        done = true;
      }
      else
      {
        // If the solver has refined, we do not allow it to double its 
	// timestep for anoter step
        if (justrefined)
	{
          justrefined = false;
        }
	else
	{
        
	  double tmp = 2.0*ldt;
          //if (fabs(ldt-tmp)<eps)
          if (tmp >= _ldt)
            ldt = _ldt;
          else
            ldt = tmp;
        }

	// If we are passed t_end
        if ((t + ldt) > t_end)
          ldt = t_end - t;
      }

      // Add increment
      for (i = 0; i < _ode->num_states(); ++i)
        y[i] += z1[i];
#ifdef DEBUG
      newton_accepted1.push_back(1);
#endif    
    }
    else
    {
      ldt /= 2.0;
      justrefined = true;
#ifdef DEBUG
      newton_accepted1.push_back(0);
#endif    
    }
  }
#ifdef DEBUG
  // Lower level than DEBUG!
  log(5, "ThetaSolver done with comp_jac = %d and rejected = %d at t=%1.2e in " \
      "%ld steps\n", jac_comp, rejects, t, num_tsteps);
#endif

}
//-----------------------------------------------------------------------------



