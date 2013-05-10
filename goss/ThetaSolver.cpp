#include <stdio.h>
#include <cmath>

#include "ThetaSolver.h"


using namespace goss;

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
  
  justrefined = false;
  num_tsteps = 0;

  stages = 1;
  newton_iter1.clear();
  newton_accepted1.clear();
  dt_v.clear();
  
  ImplicitODESolver::reset();
}
//-----------------------------------------------------------------------------
void ThetaSolver::forward(double* y, double t, double interval)
{

  uint i;
  double t_end = t + interval;
  _dt = _ldt > 0 ? _ldt : interval;
  
  for (i = 0; i < _ode->num_states(); ++i)
    _prev[i] = 0.0;

  bool step_ok, done = false;

  // A way to check if we are at t_end.
  const double eps = 1e-14;

  while (!done)
  {
    num_tsteps += 1;
  
    if (recompute_jacobian)
    {
      _ode->compute_jacobian(t, y, &jac[0]);

      // Build Theta discretization of jacobian
      mult(-_dt*(1-theta), &jac[0]);
      add_identity(&jac[0]);

      // Factorize jacobian
      _ode->lu_factorize(&jac[0]);
      jac_comp += 1;
      
    }

    // Use 0.0 z1:
    for (i = 0; i < _ode->num_states(); ++i)
      z1[i] = 0.0;

    // Explicit eval
    _ode->eval(y, t, &ft1[0]);
    
    for (i = 0; i < _ode->num_states(); ++i)
      _prev[i] = theta*ft1[i];
    
    // Solve for increment
    step_ok = newton_solve(&z1[0], &_prev[0], y, t + _dt, _dt, theta);    
#ifdef DEBUG
    newton_iter1.push_back(newtonits);
    dt_v.push_back(_dt);
#endif

    // Newton step OK
    if (step_ok)
    {
    
      t+=_dt;
      
      if (std::fabs(t-t_end)<eps)
      {
        done = true;
        //printf("Done at t=%1.4e\n",t);
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
        
	  double tmp = 2.0*_dt;
          //if (fabs(ldt-tmp)<eps)
          if (tmp>=_ldt)
            _dt = _ldt;
          else
            _dt = tmp;
        }

	// If we are passed t_end
        if ((t + _dt) > t_end)
          _dt = t_end - t;
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
      _dt /= 2.0;
      justrefined = true;
#ifdef DEBUG
      newton_accepted1.push_back(0);
#endif    
    }
  }
#ifdef DEBUG
  printf("ThetaSolver done with comp_jac = %d and rejected = %d at t=%1.2e in "\
	 "%ld steps\n", jac_comp, rejects, t, num_tsteps);
#endif

}
//-----------------------------------------------------------------------------


