#include <cassert>
#include <cmath>

#include "ImplicitEuler.h"

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
  z1(new double[solver.num_states()]), justrefined(solver.justrefined)
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
  z1.reset(new double[num_states()]);

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
  const double eps = 1e-14;

  while (!done)
  {
    num_tsteps += 1;
    
    if (recompute_jacobian)
    {
      compute_jacobian(t, y);

      // Build Euler discretization of jacobian
      mult(-_dt, jac);
      add_identity(jac);

      // Factorize jacobian
      lu_factorize(jac);
      jac_comp += 1;
    }

    // Use 0.0 as initial guess
    for (i = 0; i < num_states(); ++i)
      z1[i] = 0.0;

    // Solve for increment
    step_ok = newton_solve(z1.get(), _prev.get(), y, t + _dt, _dt, 1.0);    
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
	// timestep for anoter step
        if (justrefined)
	{
          justrefined = false;
        }
	else
	{
          double tmp = 2.0*_dt;
          //if (fabs(ldt-tmp) < eps)
          if (tmp >= _ldt)
            _dt = _ldt;
          else
            _dt = tmp;
        }
	
	// If we are passed t_end
        if ((t + _dt) > t_end)
          _dt = t_end - t;
	
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
      //std::cout << "Newton step NOT OK: " << std::endl;
      _dt /= 2.0;
      justrefined = true;
#ifdef DEBUG
      newton_accepted1.push_back(0);
#endif    
    }
  }
#ifdef DEBUG
  printf("ImplicitEuler done with comp_jac = %d and rejected = %d at t=%1.2e in "\
	 "%ld steps\n", jac_comp, rejects, t, num_tsteps);
#endif
}
//-----------------------------------------------------------------------------


