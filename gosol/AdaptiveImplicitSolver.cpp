#include "AdaptiveImplicitSolver.h"
#include <math.h>
#include <iostream>
#include <stdio.h>

using namespace gosol;

//-----------------------------------------------------------------------------
AdaptiveImplicitSolver:: AdaptiveImplicitSolver (ODE* ode_) 
  //-----------------------------------------------------------------------------
{ 
  attach(ode_); 
  //init(); 
} 

//-----------------------------------------------------------------------------
AdaptiveImplicitSolver:: AdaptiveImplicitSolver() 
  //-----------------------------------------------------------------------------
{
  //printf("AdaptiveImplicitSolver ()\n"); 
  //init(); 
}

//-----------------------------------------------------------------------------
AdaptiveImplicitSolver:: ~AdaptiveImplicitSolver ()
  //-----------------------------------------------------------------------------
{
  //printf("~AdaptiveImplicitSolver ()\n");
  //free(blog);
  //free(tvec);
};

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: init()
  //-----------------------------------------------------------------------------
{
  printf("AdaptiveImplicitSolver::init\n");
  single_step_mode=false;
  atol    = 1.0e-5;
  rtol    = 1.0e-8;
  facmin  = 0.5; // We can not choose the next timestep more then
  // half of the previous timestep
  facmaxb = 2.0; // We can not choose the next timestep more then 
  // double of the previous timestep
  facmax  = facmaxb;
  stabfac = 0.9;//pow(0.25,1/(iord+1));
  stabdown = 1.0;
  stabup   = 1.2;
  err_old  = -1.0;
  reached_tend = false;
  num_accepted=0;
  num_rejected=0;
  ImplicitODESolver::init();
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: setTol(double atol_, double rtol_)
  //-----------------------------------------------------------------------------
{
  atol=atol_;
  rtol=rtol_;
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: setIord(int iord_)
  //-----------------------------------------------------------------------------
{
  iord=iord_;
}

/*
//-----------------------------------------------------------------------------
void  AdaptiveImplicitSolver:: attach(ODE* ode_)
  //-----------------------------------------------------------------------------
{
  ode = ode_;
  ode->size() = ode->evalSize();
}
*/

//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver:: dtinit(double t, double* y0, double* y1, double* f0, double* f1, double iord)
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
  //ode->eval(y0, t, f0);
  for (i=0; i<ode->size();++i){
    sk   = atol + rtol*fabs(y0[i]);
    tmp  = f0[i]/sk;
    dnf += tmp*tmp;
    tmp  = y0[i]/sk;
    dny += tmp*tmp;
    //dnf += pow(f0[i]/sk,2);
    //dny += pow(y0[i]/sk,2); 
  }
  if (dnf <= 1.0e-10 || dny <= 1.0e-10){
    dt = 1.0e-6;
  }else{
    dt = 0.01*sqrt(dny/dnf);
  }

  // Should we have a dt_max??

  // Perform an explicit Euler step
  for (i=0; i<ode->size();++i){
    y1[i] = y0[i] + dt*f0[i];
  }
  ode->eval(y1, t+dt, f1);

  // Estimate the second derivative of the solution
  double der2 = 0.0;
  for (i=0;i<ode->size();++i){
    sk    = atol + rtol*fabs(y1[i]);
    tmp   = ((f1[i]-f0[i])/sk);
    der2 += tmp*tmp;
    //der2 += pow(((f1[i]-f0[i])/sk),2);
  }
  der2 = sqrt(der2)/dt;

  // Step size is computed such that
  // dt**iord*max(norm(f0),norm(der2)) = 0.01

  double der12;
  if (fabs(der2)>= sqrt(dnf)){
    der12=fabs(der2);
  }else{
    der12=sqrt(dnf);
  }
  double dt1;
  if (der12 <= 1.0e-15){
    if (1.0e-6 > fabs(dt)*1.0e-3){
      dt1 = 1.0e-6;
    }else{
      dt = fabs(dt)*1.0e-3;
    }
  }else{
    dt1 = pow((0.01/der12),(1.0/iord));
  }



  if (100*fabs(dt)<dt1){
    return 100*fabs(dt);
  }else{
    return dt1;
  }
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: newTimeStep(double* y, double* yn, double* e, double t_end)
//-----------------------------------------------------------------------------
{
  double err = 0.0;
  dt_prev=dt;
  double sk;
  //bool done = false;
  double eps = 1e-14;// A way to chech if we are at t_end.
  double max;
  double yi_abs, yni_abs, tmp;
  recompute_jacobian = true;
  for (int i=0; i< ode->size(); ++i) {
    yi_abs = fabs(y[i]);
    yni_abs = fabs(yn[i]);
    max = (yi_abs > yni_abs ? yi_abs : yni_abs);
    sk   = atol+rtol*max;
    tmp = e[i]/sk;
    err += tmp*tmp;
  }
  
  err=sqrt(err/ode->size());

  // If the error is smaller then 1, the timestep is accepted, and we advance
  // If not, the timestep is rejected
  if (err<=1){
    t+=dt;
    num_accepted+=1;
    step_accepted=true;
    //std::cout << "fabs(t-t_end)="<<fabs(t-t_end)<<", t, t_end = " << t << "," << t_end <<std::endl;
    if (fabs(t-t_end)<eps){
      //std::cout << "done=true" << std::endl;
      reached_tend=true;
    }
  } else {
    num_rejected += 1;
    step_accepted = false;
  }

  // Computation of dtnew
  double lstabfac=stabfac*(2*maxits+1)/((double)(2*maxits+newtonits));
  //printf("lstabfac=%1.2e\n",lstabfac);
  double fac = lstabfac*pow((1.0/err),(1.0/(iord+1)));
  if (facmin>fac){
    fac=facmin;
  }else if (fac>facmax){
    fac=facmax;
  }
  //if the timestep i rejected, we prevent the next timestep from increasing
  if (!step_accepted){
    facmax=1.0;
  }else{
    facmax=facmaxb;
  }
  if (err_old>0){
    fac*=dt/dt_old*pow((err_old/err),(1.0/(iord+1)));
  }
  if (fac<stabup && fac>stabdown){
    //printf("frac=%1.2e\n",fac);
    fac=1.0;
    recompute_jacobian = false;
  }
  dt*=fac;

  //std::cout << "t+dt="<<t+dt<<std::endl;
  ldt=dt;// Saves the timestep to be used as initial guess for next macro step
  if (t+dt>=t_end){
    dt=t_end-t;
  }

  dt_old = dt_prev;
  err_old = err;
}


#ifdef DEBUG
//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: logData(double dt, bool accepted)
  //-----------------------------------------------------------------------------
{
  dt_v.push_back(dt);
  accept_v.push_back(accepted);
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: dtVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  //std::cout << "Size = " << dt_v.size() << std::endl;
  //std::cout << "dt_v.end()="<<dt_v[dt_v.size()-2]<<std::std::endl;
  res->n    = dt_v.size();
  res->data = static_cast<double*>(malloc(sizeof(double)*(dt_v.size())));
  for( unsigned int i = 0; i < dt_v.size(); ++i ) {
    res->data[i] = dt_v[i];
  }
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: acceptedVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  res->n    = accept_v.size();
  res->data = static_cast<double*>(malloc(sizeof(double)*(accept_v.size())));
  for( unsigned int i = 0; i < accept_v.size(); ++i ) {
    res->data[i] = float(accept_v[i]);
  }
}

//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::getCurrentTime()
  //-----------------------------------------------------------------------------
{
  return t;
}

//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::getCurrentTimeStep()
  //-----------------------------------------------------------------------------
{
  return dt_prev;
}

#else
//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::getCurrentTime()
  //-----------------------------------------------------------------------------
{
  printf("NOT IN DEBUG MODE\n");
  return t;
}

//-----------------------------------------------------------------------------
double AdaptiveImplicitSolver::getCurrentTimeStep()
  //-----------------------------------------------------------------------------
{
  printf("NOT IN DEBUG MODE\n");
  return dt; 
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: logData(double dt, bool accepted)
  //-----------------------------------------------------------------------------
{
  printf("DEBUG OFF!\n");
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: dtVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  printf("DEBUG OFF!\n");
}

//-----------------------------------------------------------------------------
void AdaptiveImplicitSolver:: acceptedVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  printf("DEBUG OFF!\n");
}

#endif




