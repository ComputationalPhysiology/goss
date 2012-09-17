#include "AdaptiveImplicitSolver.h"
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace goss;

//-----------------------------------------------------------------------------
AdaptiveImplicitSolver::AdaptiveImplicitSolver (ODE* ode_) 
{ 
  attach(ode_); 
  //init(); 
} 

//-----------------------------------------------------------------------------
AdaptiveImplicitSolver::AdaptiveImplicitSolver() 
{
  //printf("AdaptiveImplicitSolver ()\n"); 
  //init(); 
}

//-----------------------------------------------------------------------------
AdaptiveImplicitSolver::~AdaptiveImplicitSolver ()
{
  //printf("~AdaptiveImplicitSolver ()\n");
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver::init()
{
  // FIXME: Flesh out constants and initialize in constructor
  printf("AdaptiveImplicitSolver::init\n");
  single_step_mode=false;
  _atol    = 1.0e-5;
  _rtol    = 1.0e-8;
  
  // We can not choose the next timestep more then half of the previous 
  // timestep
  facmin  = 0.5; 
  
  // We can not choose the next timestep more then double of the previous 
  // timestep  
  facmaxb = 2.0; 
  facmax  = facmaxb;
  stabfac = 0.9; //std::pow(0.25,1/(iord+1));
  stabdown = 1.0;
  stabup   = 1.2;
  err_old  = -1.0;
  reached_tend = false;

  num_accepted = 0;
  num_rejected = 0;
  ImplicitODESolver::init();
}
// FIXME: Have a closer look of origin of pointers passed to this function.
//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::dtinit(double t, double* y0, double* y1, 
				      double* f0_, double* f1_, double iord)
{
  // Computation of an initial step size guess
  //
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

  // FIXME: Why is it not evaluated?
  //_ode->eval(y0, t, f0_);
  for (i=0; i<ode_size(); ++i)
  {
    sk   = _atol + _rtol*std::fabs(y0[i]);
    tmp  = f0_[i]/sk;
    dnf += tmp*tmp;
    tmp  = y0[i]/sk;
    dny += tmp*tmp;
    //dnf += std::pow(f0_[i]/sk, 2);
    //dny += std::pow(y0[i]/sk, 2); 
  }
  
  if (dnf <= 1.0e-10 || dny <= 1.0e-10)
    dt = 1.0e-6;
  else
    dt = 0.01*std::sqrt(dny/dnf);

  // Should we have a dt_max??

  // Perform an explicit Euler step
  for (i = 0; i < ode_size(); ++i)
    y1[i] = y0[i] + dt*f0_[i];

  _ode->eval(y1, t+dt, f1_);

  // Estimate the second derivative of the solution
  double der2 = 0.0;
  for (i = 0; i < ode_size(); ++i)
  {
    sk    = _atol + _rtol*std::fabs(y1[i]);
    tmp   = ((f1_[i]-f0_[i])/sk);
    der2 += tmp*tmp;
    //der2 += std::pow(((f1_[i]-f0_[i])/sk), 2);
  }
  der2 = std::sqrt(der2)/dt;

  // Step size is computed such that
  // dt**iord*max(norm(f0_),norm(der2)) = 0.01

  double der12;
  if (std::fabs(der2) >= std::sqrt(dnf))
    der12 = std::fabs(der2);
  else
    der12 = std::sqrt(dnf);
  
  double dt1;
  if (der12 <= 1.0e-15)
  {
    if (1.0e-6 > std::fabs(dt)*1.0e-3)
      dt1 = 1.0e-6;
    else
      dt1 = std::fabs(dt)*1.0e-3;
  }
  else
  {
    dt1 = std::pow((0.01/der12), (1.0/iord));
  }

  if (100*std::fabs(dt) < dt1)
    return 100*std::fabs(dt);
  
  return dt1;
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver::new_time_step(double* y, double* yn, double* e, double t_end)
{
  uint i;
  double err = 0.0;
  _dt_prev = _dt;
  //bool done = false;

  // A way to check if we are at t_end.
  const double eps = 1e-14;
  recompute_jacobian = true;

  double yi_abs, yni_abs, max, sk, tmp;
  for (i = 0; i < ode_size(); ++i) 
  {
    yi_abs = std::fabs(y[i]);
    yni_abs = std::fabs(yn[i]);
    max = std::max(yi_abs, yni_abs);//(yi_abs > yni_abs ? yi_abs : yni_abs);
    sk = _atol + _rtol*max;
    tmp = e[i]/sk;
    err += tmp*tmp;
  }
  
  err = std::sqrt(err/ode_size());

  // If the error is smaller than 1, the timestep is accepted, and we advance
  // If not, the timestep is rejected
  if (err <= 1)
  {
    _t += _dt;
    num_accepted += 1;
    step_accepted = true;
    //std::cout << "std::fabs(t-t_end)="<<std::fabs(t-t_end)<<", t, t_end = " << t << "," << t_end <<std::endl;
    if (std::fabs(_t - t_end) < eps)
    {
      //std::cout << "done=true" << std::endl;
      reached_tend = true;
    }
  } 
  else 
  {
    num_rejected += 1;
    step_accepted = false;
  }

  // Computation of dtnew
  const double lstabfac = stabfac*(2*maxits+1)/((2.0*maxits+newtonits));
  
  //printf("lstabfac=%1.2e\n", lstabfac);
  double fac = lstabfac*std::pow((1.0/err), (1.0/(_iord+1)));
  if (facmin > fac)
    fac = facmin;

  else if (fac > facmax)
    fac = facmax;

  //if the timestep i rejected, we prevent the next timestep from increasing
  if (!step_accepted)
    facmax = 1.0;
  else
    facmax = facmaxb;

  if (err_old > 0)
    fac *= _dt/dt_old*std::pow((err_old/err), (1.0/(_iord + 1)));
  
  if (fac < stabup && fac > stabdown){
    //printf("frac=%1.2e\n",fac);
    fac = 1.0;
    recompute_jacobian = false;
  }
  _dt *= fac;

  //std::cout << "t+dt="<<t+dt<<std::endl;
  _ldt = _dt;// Saves the timestep to be used as initial guess for next macro step

  if (_t + _dt >= t_end)
    _dt = t_end - _t;
  
  dt_old = _dt_prev;
  err_old = err;
}
//-----------------------------------------------------------------------------
#ifdef DEBUG
void AdaptiveImplicitSolver::log_data(double dt, bool accepted)
{
  dt_v.push_back(dt);
  accept_v.push_back(accepted);
}
//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver::dt_vector(DoubleVector *res)
{
  res->n = dt_v.size();
  res->data = new double[dt_v.size()]; 
  for(uint i = 0; i < dt_v.size(); ++i)
    res->data[i] = dt_v[i];
}
//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver::accepted_vector(DoubleVector *res)
{
  res->n    = accept_v.size();
  res->data = new double[accept_v.size()]; 
  for(uint i = 0; i < accept_v.size(); ++i)
    res->data[i] = float(accept_v[i]);
}
//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::get_current_time()
{
  return _t;
}
//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::get_current_time_step()
{
  return _dt_prev;
}
#else
//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver::log_data(double, bool)
{
  std::cout << "DEBUG OFF!" << std::endl; 
}
//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver::dt_vector(DoubleVector*)
{
  std::cout << "DEBUG OFF!" << std::endl; 
}
//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver::accepted_vector(DoubleVector*)
{
  std::cout << "DEBUG OFF!" << std::endl; 
}
//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::get_current_time()
{
  printf("NOT IN DEBUG MODE\n");
  return _t;
}
//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::get_current_time_step()
{
  printf("NOT IN DEBUG MODE\n");
  return _dt_prev;
}
//-----------------------------------------------------------------------------
#endif




