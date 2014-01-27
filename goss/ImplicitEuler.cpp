#include <cassert>
#include <cmath>
#include <cstdio>

#include "ImplicitEuler.h"
#include "log.h"
#include "constants.h"

using namespace goss;

//-----------------------------------------------------------------------------
ImplicitEuler::ImplicitEuler() : ImplicitODESolver(), newton_iter1(0), 
				 newton_accepted1(0), dt_v(0), z1(0), 
				 justrefined(false)
{ 
  // Do nothing
} 
//-----------------------------------------------------------------------------
ImplicitEuler::ImplicitEuler(boost::shared_ptr<ODE> ode, double ldt) : 
  ImplicitODESolver(ldt), newton_iter1(0), newton_accepted1(0), dt_v(0), 
  z1(0), justrefined(false)
{ 
  attach(ode);
} 
//-----------------------------------------------------------------------------
ImplicitEuler::ImplicitEuler(double ldt) : ImplicitODESolver(ldt), newton_iter1(0), 
					   newton_accepted1(0), dt_v(0), z1(0), 
					   justrefined(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitEuler::ImplicitEuler(const ImplicitEuler& solver) : 
  ImplicitODESolver(solver), newton_iter1(solver.newton_iter1), 
  newton_accepted1(solver.newton_accepted1), dt_v(solver.dt_v), 
  z1(solver.num_states()), justrefined(solver.justrefined)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitEuler::~ImplicitEuler()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ImplicitEuler::attach(boost::shared_ptr<ODE> ode)
{

  // Attach ode using bases
  ImplicitODESolver::attach(ode);

  // Init memory
  z1.resize(num_states());

}
//-----------------------------------------------------------------------------
void ImplicitEuler::reset()
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
void ImplicitEuler::compute_factorized_jacobian(double* y, double t, double dt)
{
  
  // Let ODE compute the jacobian
  _ode->compute_jacobian(y, t, &jac[0]);

  // Build Euler discretization of jacobian
  mult(-dt, &jac[0]);
  add_mass_matrix(&jac[0]);

  // Factorize the jacobian
  _ode->lu_factorize(&jac[0]);
  jac_comp += 1;

}
//-----------------------------------------------------------------------------
void ImplicitEuler::forward(double* y, double t, double interval) 
{

  assert(_ode);

  uint i;
  const double t_end = t + interval;

  _dt = _ldt > 0 ? _ldt : interval;
  
  for (i = 0; i < num_states(); ++i)
    _prev[i] = 0.0;

  bool step_ok, done = false;

  // A way to check if we are at t_end.
  const double eps = GOSS_EPS*1000;

  while (!done)
  {
  
    num_tsteps += 1;
    
    // Recompute the jacobian if nessesary
    if (recompute_jacobian)
    {
      compute_factorized_jacobian(y, t, _dt);
    }

    // Use 0.0 as initial guess
    for (i = 0; i < num_states(); ++i)
      z1[i] = 0.0;

    // Solve for increment
    step_ok = newton_solve(&z1[0], &_prev[0], y, t + _dt, _dt, 1.0);    
#ifdef DEBUG
    newton_iter1.push_back(newtonits);
    dt_v.push_back(_dt);
#endif
    
    // Newton step OK
    if (step_ok)
    {
      //std::cout << "Newton step OK: " << std::endl;
      t += _dt;
      if (std::fabs(t - t_end) < eps)
      {
        done = true;
        //printf("Done at t=%1.4e\n",t);
      }
      else
      {
        // If the solver has refined, we do not allow it to double its 
	// timestep for another step
        if (justrefined)
	{
          justrefined = false;
        }
	else
	{
          const double tmp = 2.0*_dt;
	  //if (fabs(_ldt-tmp) < eps)
	  if (_ldt > 0 && tmp >= _ldt)
            _dt = _ldt;
          else
            _dt = tmp;
	  //if (std::abs(_dt-tmp/2)<GOSS_EPS)
	    goss_debug2("Changing dt from %e to %e", tmp/2, _dt);
        }
	
	// If we are passed t_end
        if ((t + _dt) > t_end)
	{
	  _dt = t_end - t;
	  goss_debug1("Adapting timestep due to t_end: dt %e", _dt);
	}
      }
      
      // Add increment
      for (i = 0; i < num_states(); ++i)
        y[i] += z1[i];
#ifdef DEBUG
      newton_accepted1.push_back(1);
#endif    
    }
    else
    {
      _dt /= 2.0;
      goss_debug1("Reducing dt: %e", _dt);

      recompute_jacobian = true;
      justrefined = true;
#ifdef DEBUG
      newton_accepted1.push_back(0);
#endif    
    }
  }
#ifdef DEBUG
  // Lower level than DEBUG!
  log(5, "ImplicitEuler done with comp_jac = %d and rejected = %d at t=%1.2e in " \
      "%ld steps\n", jac_comp, rejects, t, num_tsteps);
#endif
}
//-----------------------------------------------------------------------------


