#include "RKF32.h"
#include "math.h"
#include <iostream>
#include <cstring>
#include <cstdlib>

using namespace goss;

//-----------------------------------------------------------------------------
RKF32:: RKF32(ODE* ode, double ldt) 
//-----------------------------------------------------------------------------
  : AdaptiveExplicitSolver(ode, ldt, 0.0), 
    ki(0), k1(0), k2(0), k3(0), k4(0), yn(0), e(0)
{ 
  init();
  attach(ode);
}

//-----------------------------------------------------------------------------
RKF32:: RKF32(double ldt) 
//-----------------------------------------------------------------------------
  : AdaptiveExplicitSolver(ldt), ki(0), k1(0), k2(0), k3(0), k4(0), yn(0), e(0)
{
  init();
  AdaptiveExplicitSolver::init();
}

//-----------------------------------------------------------------------------
void RKF32:: init() 
  //-----------------------------------------------------------------------------
{
  iord = 3;
}

//-----------------------------------------------------------------------------
RKF32:: ~RKF32() 
  //-----------------------------------------------------------------------------
{
  swap = NULL;
  free(ki);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(yn);
  free(e);
}

//-----------------------------------------------------------------------------
void  RKF32:: attach(ODE* ode_)
  //-----------------------------------------------------------------------------
{
  swap = NULL;
  free(ki);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(yn);
  free(e);

  ode = ode_;
  ki   = static_cast<double*>(malloc(ode->size()*sizeof*ki));
  k1   = static_cast<double*>(malloc(ode->size()*sizeof*k1));
  k2   = static_cast<double*>(malloc(ode->size()*sizeof*k2));
  k3   = static_cast<double*>(malloc(ode->size()*sizeof*k3));
  k4   = static_cast<double*>(malloc(ode->size()*sizeof*k4));
  yn   = static_cast<double*>(malloc(ode->size()*sizeof*yn));
  e    = static_cast<double*>(malloc(ode->size()*sizeof*e));
  swap = 0;

  first = true;

  a21  = 1.0/2.0;
  a32  = 3.0/4.0;
  b1   = 2.0/9.0;
  b2   = 1.0/3.0;
  b3   = 4.0/9.0;
  bh1  = 7.0/24.0;
  bh2  = 1.0/4.0;
  bh3  = 1.0/3.0; 
  bh4  = 1.0/8.0;
  d1   = b1-bh1;
  d2   = b2-bh2;
  d3   = b3-bh3;
  d4   =   -bh4;
  c2   = 1.0/2.0;
  c3   = 3.0/4.0;
  nbytes = ode->size()*sizeof(double);
  nfevals =0;
  ndtsa   =0;
  ndtsr   =0;
}


//-----------------------------------------------------------------------------
void RKF32:: forward(double* y, double t_, double interval) {
  //-----------------------------------------------------------------------------
  t = t_;
  int i;
  double t_end;

  t_end = t + interval;

  reached_tend = false;

  // We swap the result vektor if a timestep is accepted. We therefore need to store the pointer to the initial y-vector in order to ensure that this memory segment contains the final result when the end is reached ...
  retPtr = y;

  if (first) {
    ode->eval(y, t, k1);
    nfevals += 1;
  }

  if (ldt < 0) {
    dt = dtinit(t, y, yn, k1, k2, iord);
    nfevals += 1;
  } else {
    dt = ldt;
  }
#ifdef DEBUG
  // log data
  dt_v.push_back(dt);
#endif

  while (!reached_tend){
    for (i = 0; i < ode->size(); ++i)
      ki[i] = y[i] + dt*a21*k1[i];

    ode->eval(ki, t + c2*dt, k2);
    for (i = 0; i < ode->size(); ++i)
      ki[i] = y[i] + dt*a32*k2[i];

    ode->eval(ki, t + c3*dt, k3);
    // We assemble the new y
    for (i = 0; i < ode->size(); ++i)
      yn[i] = y[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i]);
    //yn[i] = y[i] + dt*(bh1*k1[i]+bh2*k2[i]+bh3*k3[i]+bh4*k4[i]);

    // We compute the first quadrature node for the next iteration (FSAL)
    ode->eval(yn, t+dt, k4);
    nfevals += 3;
    // We compute the error vector
    for (i=0; i < ode->size(); ++i)
      e[i] = dt*(d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i]);

    newTimeStep(y, yn, e, t_end);

#ifdef DEBUG
    logData(dt,step_accepted);
#endif
    if (step_accepted){
      ndtsa += 1;
      swap = y;
      y    = yn;
      yn   = swap;
      swap = k1;
      k1   = k4;
      k4   = swap;
#ifdef DEBUG
      if (single_step_mode){
        if (retPtr!=y){
          memcpy(retPtr, y, nbytes);
          yn = y;
        }
        swap = 0;
        return;
      }
#endif
    } else
      ndtsr += 1;

  }

  // This is a copy to ensure that the input Ptr contains the final solution
  // This can probably be done in a more elegant way
  if (retPtr!=y){
    memcpy(retPtr, y, nbytes);
    yn = y;
  }
#ifdef DEBUG
  dt_v.pop_back();
#endif
}

#ifdef DEBUG
//-----------------------------------------------------------------------------
void RKF32:: logData(double dt, bool accepted)
  //-----------------------------------------------------------------------------
{
  dt_v.push_back(dt);
  accept_v.push_back(accepted);
}

//-----------------------------------------------------------------------------
void RKF32:: dtVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  res->n    = dt_v.size();
  res->data = static_cast<double*>(malloc(sizeof(double)*(dt_v.size())));
  for( unsigned int i = 0; i < dt_v.size(); ++i ) {
    res->data[i] = dt_v[i];
  }
}

//-----------------------------------------------------------------------------
void RKF32:: acceptedVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  res->n    = accept_v.size();
  res->data = static_cast<double*>(malloc(sizeof(double)*(accept_v.size())));
  for( unsigned int i = 0; i < accept_v.size(); ++i ) {
    res->data[i] = float(accept_v[i]);
  }
}

#else
//-----------------------------------------------------------------------------
void RKF32:: logData(double dt, bool accepted)
  //-----------------------------------------------------------------------------
{
  std::cout << "DEBUG OFF!" << std::endl; 
}

//-----------------------------------------------------------------------------
void RKF32:: dtVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  std::cout << "DEBUG OFF!" << std::endl; 
}

//-----------------------------------------------------------------------------
void RKF32:: acceptedVector(DoubleVector *res)
  //-----------------------------------------------------------------------------
{
  std::cout << "DEBUG OFF!" << std::endl; 
}

#endif
