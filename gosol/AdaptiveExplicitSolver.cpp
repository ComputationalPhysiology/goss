#include "AdaptiveExplicitSolver.h"
#include <math.h>
#include <iostream>
#include <stdio.h>

using namespace gosol;

//-----------------------------------------------------------------------------
AdaptiveExplicitSolver:: AdaptiveExplicitSolver (ODE* ode, double ldt, double dt)
//-----------------------------------------------------------------------------
  : ODESolver(ode, ldt, dt), dt_v(0), accept_v(0)
{ 
  init();
} 

//-----------------------------------------------------------------------------
AdaptiveExplicitSolver:: AdaptiveExplicitSolver(double ldt, double dt) 
//-----------------------------------------------------------------------------
  : ODESolver(ldt, dt), dt_v(0), accept_v(0)
{
  init();
}


//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver:: init()
  //-----------------------------------------------------------------------------
{
  single_step_mode = false;
  reached_tend     = false;
  num_accepted = 0;
  num_rejected = 0;
  atol    = 1.0e-5;
  rtol    = 1.0e-8;
  itol    = 0;
  facmin  = 0.5; // We can not choose the next timestep more then
  // half of the previous timestep
  facmaxb = 2.0; // We can not choose the next timestep more then 
  // double of the previous timestep
  facmax  = facmaxb;
  stabfac = 0.9;//pow(0.25,1/(iord+1));
}

//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver:: setTol(double atol_, double rtol_)
  //-----------------------------------------------------------------------------
{
  atol = atol_;
  rtol = rtol_;
}

//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver:: setIord(int iord_)
  //-----------------------------------------------------------------------------
{
  iord = iord_;
}


//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver:: dtinit(double t, double* y0, double* y1, double* f0, double* f1, double iord)
  //-----------------------------------------------------------------------------
{
  // Computation of an initial step size guess
  //
  // Compute a first guess for explicit Euler as
  // H = 0.01 * norm(y0)/(norm(f0)
  // The increment for explicit Euler is small compared to the solution
  // We assume that y0 and f0 are computed.
  // y1 and f1 are just pointers to contigous memory which this 
  // function borrows

  int i;
  double dnf = 0.0;
  double dny = 0.0;
  double sk,dt,tmp;
  if (itol == 0){
    for (i=0; i<ode->size();++i){
      sk   = atol + rtol*fabs(y0[i]);
      tmp  = f0[i]/sk;
      dnf += tmp*tmp;
      tmp  = y0[i]/sk;
      dny += tmp*tmp;
    }
  } else {
    std::cout << " Not implemented yeat " << std::endl;
    //for (i=0; i<ode->size();++i){
    //    sk   = atol[i] + rtol[i]*fabs(y0[i]);
    //    dnf += pow(f0[i]/sk,2);
    //    dny += pow(y0[i]/sk,2); 
    //}
    //for i in xrange(0,n):
    //    sk   = atol[i] + rtol[i]*math.fabs(y0[i])
    //    dnf += (f0[i]/sk)**2
    //    dny += (y0[i]/sk)**2
    //print "dnf = ", dnf, " dny = ", dny
  }
  if (dnf <= 1.0e-10 || dny <= 1.0e-10)
  {
    dt = 1.0e-6;
  } 
  else
  {
    dt = 0.01*sqrt(dny/dnf);
  }

  // Should we have a dt_max??
  // Perform an explicit Euler step
  for (i=0; i<ode->size();++i){
    y1[i] = y0[i] + dt*f0[i];
  }
  ode->eval(y1, t + dt, f1);

  // Estimate the second derivative of the solution
  double der2 = 0.0;
  if (itol == 0)
  {
    for (i=0;i<ode->size();++i)
    {
      sk    = atol + rtol*fabs(y1[i]);
      tmp   = (f1[i] - f0[i])/sk;
      der2 += tmp*tmp;
    }
  }
  else
  {
    std::cout << " Not implemented yeat" << std::endl;
    //for (i=0;i<ode->size();++i){
    //    sk    = atol[i] + rtol[i]*fabs(y1[i]);
    //    der2 += pow(((f1[i]-f0[i])/sk),2);
    //}
    //for i in xrange(0,n):
    //     sk    = atol[i] + rtol[i]*math.fabs(y0[i])
    //     der2 += ((f1[i]-f0[i])/sk)**2
  }
  der2 = sqrt(der2)/dt;

  // Step size is computed such that
  // dt**iord*max(norm(f0),norm(der2)) = 0.01

  double der12;
  if (fabs(der2) >= sqrt(dnf))
  {
    der12 = fabs(der2);
  }
  else
  {
    der12 = sqrt(dnf);
  }
  double dt1;
  if (der12 <= 1.0e-15)
  {
    if (1.0e-6 > fabs(dt)*1.0e-3)
    {
      dt1 = 1.0e-6;
    } 
    else 
    {
      dt = fabs(dt)*1.0e-3;
    }
  } 
  else 
  {
    dt1 = pow((0.01/der12),(1.0/iord));
  }

  if (100*fabs(dt)<dt1)
  {
    return 100*fabs(dt);
  }
  else
  {
    return dt1;
  }
}


//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver:: newTimeStep(double* y, double* yn, double* e, double t_end)
  //-----------------------------------------------------------------------------
{
  double err = 0.0;
  double sk;
  double eps = 1e-14;// A way to chech if we are at t_end.
  double max;
  double yi_abs, yni_abs, tmp;
  dt_prev = dt;
  for (int i=0; i< ode->size(); ++i) {
    yi_abs  = fabs(y[i]);
    yni_abs = fabs(yn[i]);
    max     = (yi_abs > yni_abs ? yi_abs : yni_abs);
    sk      = atol + rtol*max;
    tmp     = e[i]/sk;
    err += tmp*tmp;
  }
  err = sqrt(err/ode->size());

  // Sanity check to ensure algorithm does not break down.
  // If we encounter strange values in the error, set it large to enforce a
  // small time step for the next iteration.
  if (isnan(err) || isinf(err))
    err = 2000.0;

  // If the error is smaller then 1, the timestep is accepted, and we advance
  // If not, the timestep is rejected
  if (err <= 1.0)
  {
    t += dt;
    num_accepted += 1;
    step_accepted = true;
    if (fabs(t-t_end)<eps)
    {
      reached_tend = true;
    }
  }
  else
  {
    num_rejected += 1;
    step_accepted = false;
  }

  // Computation of dtnew
  double fac = stabfac*pow(1.0/err, 1.0/(iord+1));
  if (facmin > fac)
  {
    fac = facmin;
  }
  if (fac > facmax)
  {
    fac = facmax;
  }

  //if the timestep is rejected, we prevent the next timestep from increasing
  if (!step_accepted)
  {
    facmax = 1.0;
  }
  else
  {
    facmax = facmaxb;
  }

  dt *= fac;
  //If it is very likely that we will reach the end with two timesteps, we set the first timestp to half the distance to avoid the last timestep been very small
  if (t+dt>=t_end)
  {
    dt = t_end-t;
  }
  else if (t+1.5*dt>=t_end)
  {
    dt=(t_end-t)/2.0;
  }
  else
  {
    ldt = dt;// Saves the timestep to be used as initial guess for next macro step
  }
  //std::cout << "dt = " << dt << std::endl;
}

//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::getCurrentTime()
  //-----------------------------------------------------------------------------
{
#ifndef DEBUG
  printf("NOT IN DEBUG MODE\n");
#endif
  return t;
}

//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::getCurrentTimeStep()
  //-----------------------------------------------------------------------------
{
#ifndef DEBUG
  printf("NOT IN DEBUG MODE\n");
  return dt; 
#else
  return dt_prev;
#endif
}

