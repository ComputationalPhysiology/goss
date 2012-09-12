#include "RKF34.h"
#include "math.h"
#include <iostream>
#include <cstring>

using namespace gossol;

void  RKF34:: attach(ODE* ode_)
{
  ode = ode_;
  ki   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  f1   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  f2f5 = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  f3   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  f4   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  yn   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  e    = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  swap = NULL;


  a21  = 1.0/4.0;
  a31  = 4.0/81.0;
  a32  = 32.0/81.0;
  a41  = 57.0/98.0;
  a42  = -432.0/343.0;
  a43  = 1053.0/686.0;
  a51  = 1.0/6.0;
  a53  = 27.0/52.0;
  a54  = 49.0/156.0;
  b1   = 1.0/6.0;
  b3   = 27.0/52.0;
  b4   = 49.0/156.0;
  bh1  = 43.0/288.0;
  bh3  = 243.0/416.0;
  bh4  = 343.0/1872.0; 
  bh5  = 1.0/12.0;
  d1   = b1-bh1;
  d3   = b3-bh3;
  d4   = b4-bh4;
  d5   =   -bh5;

  c2   = 1.0/4.0;
  c3   = 4.0/9.0;
  c4   = 6.0/7.0;
  nbytes = ode->size()*sizeof(double);
  nfevals =0;
  ndtsa   =0;
  ndtsr   =0;
}



void RKF34:: forward(double* y, double t_, double interval) {
  t = t_;
  int i;
  double t_end;

  reached_tend = false;

  retPtr = y;
  t_end = t + interval;
  ode->eval(y, t, f1);
  nfevals=1;
  if (ldt<0){
    dt = dtinit(t, y, yn, f1, f2f5, iord);
    nfevals += 1;
  } else {
    dt = ldt;
  }

#ifdef DEBUG
  // log data
  dt_v.push_back(dt);
#endif

  while (!reached_tend){
    ode->eval(y          , t       , f1);
    for (i=0; i < ode->size(); ++i)
      ki[i] = y[i] + dt*a21*f1[i];

    ode->eval(ki, t+c2*dt, f2f5);
    for (i=0; i < ode->size(); ++i)
      ki[i] = y[i] + dt*(a31*f1[i] + a32*f2f5[i]);

    ode->eval(ki, t+c3*dt, f3);
    for (i=0; i < ode->size(); ++i)
      ki[i] = y[i] + dt*(a41*f1[i] + a42*f2f5[i] + a43*f3[i]);

    ode->eval(ki, t+c4*dt, f4);
    for (i=0; i < ode->size(); ++i)
      ki[i] = y[i] + dt*(a51*f1[i] + a53*f3[i] + a54*f4[i]);

    ode->eval(ki, t+c4*dt, f2f5);
    nfevals += 5;


    // We assemble the new y
    for (i=0; i < ode->size(); ++i)
      //yn[i] = y[i] + dt*(b1*f1[i]+b3*f3[i]+b4*f4[i]);// This is the third order solution
      yn[i] = y[i] + dt*(bh1*f1[i] + bh3*f3[i] + bh4*f4[i] + bh5*f2f5[i]);// This is the fourth order solution


    // We assemble the error vector
    for (i=0; i < ode->size(); ++i)
      e[i] =  dt*(d1*f1[i] + d3*f3[i] + d4*f4[i] + d5*f2f5[i]);

    newTimeStep(y, yn, e, t_end);


#ifdef DEBUG
    logData(dt,step_accepted);
#endif
    if (step_accepted){
      swap = y;
      y    = yn;
      yn   = swap;
      ndtsa += 1;
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
    }else
      ndtsr += 1;
  }
  if (retPtr!=y){
    memcpy(retPtr, y, nbytes);
    yn = y;
    swap = 0;
  }
#ifdef DEBUG
  // The last timestep suggested is 0, because we are at the end. This should therefore not be logged.
  dt_v.pop_back();
  if (single_step_mode)
    return;
#endif
}


#ifdef DEBUG
void RKF34:: logData(double dt, bool accepted) {
  dt_v.push_back(dt);
  accept_v.push_back(accepted);
}

void RKF34:: dtVector(DoubleVector *res){
  res->n    = dt_v.size();
  res->data = static_cast<double*>(malloc(sizeof(double)*(dt_v.size())));
  for( unsigned int i = 0; i < dt_v.size(); ++i ) {
    res->data[i] = dt_v[i];
  }
}

void RKF34:: acceptedVector(DoubleVector *res){
  res->n    = accept_v.size();
  res->data = static_cast<double*>(malloc(sizeof(double)*(accept_v.size())));
  for( unsigned int i = 0; i < accept_v.size(); ++i ) {
    res->data[i] = float(accept_v[i]);
  }
}
#else
void RKF34:: logData(double dt, bool accepted) {
  std::cout << "DEBUG OFF!" << std::endl; 
}

void RKF34:: dtVector(DoubleVector *res){
  std::cout << "DEBUG OFF!" << std::endl;   //std::cout << "Size = " << dt_v.size() << std::endl;
}

void RKF34:: acceptedVector(DoubleVector *res){
  std::cout << "DEBUG OFF!" << std::endl; 
}

#endif
