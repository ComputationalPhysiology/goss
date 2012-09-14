#include "AdaptiveExplicitSolver.h"
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace goss;

//-----------------------------------------------------------------------------
AdaptiveExplicitSolver::AdaptiveExplicitSolver() 
  : ODESolver(), dt_v(0), accept_v(0)
{
  init();
}
//-----------------------------------------------------------------------------
AdaptiveExplicitSolver::AdaptiveExplicitSolver (ODE* ode, double ldt, double dt)
  : ODESolver(ode, ldt, dt), dt_v(0), accept_v(0)
{ 
  init();
} 
//-----------------------------------------------------------------------------
AdaptiveExplicitSolver::AdaptiveExplicitSolver(double ldt, double dt) 
  : ODESolver(ldt, dt), dt_v(0), accept_v(0)
{
  init();
}
//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver::init()
{
  single_step_mode = false;
  reached_tend     = false;
  num_accepted = 0;
  num_rejected = 0;
  _atol    = 1.0e-5;
  _rtol    = 1.0e-8;
  _itol    = 0;

  // We can not choose the next timestep more then
  // half of the previous timestep
  facmin  = 0.5;

  // We can not choose the next timestep more then 
  // double of the previous timestep
  facmaxb = 2.0; 
  facmax  = facmaxb;
  stabfac = 0.9; // std::pow(0.25,1/(_iord+1));
}
//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::dtinit(double t, double* y0, double* y1, 
				      double* f0, double* f1, double iord)
{
  // Compute a first guess for explicit Euler as
  // H = 0.01 * norm(y0)/(norm(f0)
  // The increment for explicit Euler is small compared to the solution
  // We assume that y0 and f0 are computed.
  // y1 and f1 are just pointers to contigous memory which this 
  // function borrows

  uint i;
  double dnf = 0.0;
  double dny = 0.0;
  double sk, dt, tmp;
  if (_itol == 0)
  {
    for (i = 0; i < _ode->size(); ++i)
    {
      sk   = _atol + _rtol*fabs(y0[i]);
      tmp  = f0[i]/sk;
      dnf += tmp*tmp;
      tmp  = y0[i]/sk;
      dny += tmp*tmp;
    }
  } 
  else 
  {
    std::cout << " Not implemented yet " << std::endl;
    //for (i=0; i<_ode->size();++i){
    //    sk   = _atol[i] + _rtol[i]*fabs(y0[i]);
    //    dnf += pow(f0[i]/sk,2);
    //    dny += pow(y0[i]/sk,2); 
    //}
    //for i in xrange(0,n):
    //    sk   = _atol[i] + _rtol[i]*math.fabs(y0[i])
    //    dnf += (f0[i]/sk)**2
    //    dny += (y0[i]/sk)**2
    //print "dnf = ", dnf, " dny = ", dny
  }

  if (dnf <= 1.0e-10 || dny <= 1.0e-10)
    dt = 1.0e-6;
  else
    dt = 0.01*std::sqrt(dny/dnf);

  // Should we have a dt_max??
  // Perform an explicit Euler step
  for (i = 0; i < _ode->size(); ++i)
    y1[i] = y0[i] + dt*f0[i];

  _ode->eval(y1, t + dt, f1);

  // Estimate the second derivative of the solution
  double der2 = 0.0;
  if (_itol == 0)
  {
    for (i = 0; i < _ode->size(); ++i)
    {
      sk    = _atol + _rtol*fabs(y1[i]);
      tmp   = (f1[i] - f0[i])/sk;
      der2 += tmp*tmp;
    }
  }
  else
  {
    std::cout << " Not implemented yeat" << std::endl;
    //for (i=0;i<_ode->size();++i){
    //    sk    = _atol[i] + _rtol[i]*fabs(y1[i]);
    //    der2 += pow(((f1[i]-f0[i])/sk),2);
    //}
    //for i in xrange(0,n):
    //     sk    = _atol[i] + _rtol[i]*math.fabs(y0[i])
    //     der2 += ((f1[i]-f0[i])/sk)**2
  }
  der2 = std::sqrt(der2)/dt;

  // Step size is computed such that
  // dt**iord*max(norm(f0),norm(der2)) = 0.01

  double der12;
  if (fabs(der2) >= std::sqrt(dnf))
    der12 = fabs(der2);
  else
    der12 = std::sqrt(dnf);

  double dt1;
  if (der12 <= 1.0e-15)
  {
    if (fabs(dt)*1.0e-3 < 1.0e-6)
      dt1 = 1.0e-6;
    else 
      dt1 = fabs(dt)*1.0e-3;
  } 
  else 
    dt1 = std::pow((0.01/der12),(1.0/iord));

  if (100*fabs(dt)<dt1)
    return 100*fabs(dt);
  else
    return dt1;
}
//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver::new_time_step(double* y, double* yn, double* e, double t_end)
{
  // FIXME: Force this method to return the new time step
  double err = 0.0;
  double sk;

  // A way to check if we are at t_end.
  double eps = 1e-14;
  double max;
  double yi_abs, yni_abs, tmp;
  _dt_prev = _dt;

  for (uint i = 0; i < _ode->size(); ++i) 
  {
    yi_abs  = fabs(y[i]);
    yni_abs = fabs(yn[i]);
    max     = yi_abs > yni_abs ? yi_abs : yni_abs;
    sk      = _atol + _rtol*max;
    tmp     = e[i]/sk;
    err += tmp*tmp;
  }
  
  err = std::sqrt(err/_ode->size());

  // Sanity check to ensure algorithm does not break down.
  // If we encounter strange values in the error, set it large to enforce a
  // small time step for the next iteration.
  if (std::isnan(err) || std::isinf(err))
    err = 2000.0;

  // If the error is smaller than 1 the timestep is accepted, and we advance
  // If not, the timestep is rejected
  if (err <= 1.0)
  {
    _t += _dt;
    num_accepted += 1;
    step_accepted = true;

    if (fabs(_t - t_end) < eps)
      reached_tend = true;
  }
  else
  {
    num_rejected += 1;
    step_accepted = false;
  }

  // Computation of dtnew
  double fac = stabfac*std::pow(1.0/err, 1.0/(_iord + 1));
  if (facmin > fac)
    fac = facmin;

  if (fac > facmax)
    fac = facmax;

  // If the timestep is rejected, we prevent the next timestep from increasing
  if (!step_accepted)
    facmax = 1.0;
  else
    facmax = facmaxb;

  _dt *= fac;

  // If it is very likely that we will reach the end with two timesteps, we 
  // set the first timestp to half the distance to avoid the last timestep 
  // been very small
  if (_t + _dt >= t_end)
    _dt = t_end - _t;

  else if (_t + 1.5*_dt >= t_end)
    _dt = (t_end - _t)/2.0;

  else
    _ldt = _dt;// Saves the timestep to be used as initial guess for next macro step

  //std::cout << "dt = " << _dt << std::endl;
}
//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::get_current_time()
{
#ifndef DEBUG
  printf("NOT IN DEBUG MODE\n");
#endif
  return _t;
}

//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::get_current_time_step()
{
#ifndef DEBUG
  printf("NOT IN DEBUG MODE\n");
  return _dt; 
#else
  return _dt_prev;
#endif
}

