#include "SDIRK2O3.h"
#include <stdio.h>

using namespace gosol;

void  SDIRK2O3:: attach(ODE* ode_){
  ode = ode_;
  //    ode->attach(this);

  init();

  gamma = (3+sqrt(3))/6;
  stages = 2;
  //newton_tol = 1e-5;
  a11 = gamma;
  a21 = 1-2*gamma;
  a22 = gamma;
  b1  = 0.5;
  b2  = 0.5;
  c1  = gamma;
  c2  = 1-gamma;
  d1  = (b1-(a21*b2)/a22)/a11;
  d2  = b2/a22;
  justrefined=false;
  num_tsteps=0;


  z1 = static_cast<double*>(malloc(sizeof(double)*ode->size())); 
  z2 = static_cast<double*>(malloc(sizeof(double)*ode->size())); 
}


void SDIRK2O3:: forward(double* y, double t, double interval) {

  int i;
  double t_end = t+interval;

  dt=ldt;    

  bool step_ok, done = false;
  double eps = 1e-14;// A way to chech if we are at t_end.


  //for (j = 0; j<n_steps; ++j){
  while (!done){
    //Jacobian recomputed once for each local step
    //printf("Start of time loop\n");
    num_tsteps+=1;
    if(recompute_jacobian){
      computeJacobian(t,y);
      mult(-dt*a11,jac);
      addIdentity(jac);
      factLU(jac);
      jac_comp+=1;
#ifdef DEBUG
      printf("Recompute jac at t=%1.2e, dt = %1.4e\n",t,dt);
#endif
    }
    for (i = 0; i < ode->size(); ++i)
      prev[i] = 0.0;

    //Use y_n as initial guess for z1:
    for (i = 0; i <ode->size(); ++i)
      z1[i] = 0.0;//y[i];
    step_ok = NewtonSolve(z1,prev,y,t+c1*dt,dt,a11);    
#ifdef DEBUG
    newtonIter1.push_back(newtonits);
    dt_v.push_back(dt);
#endif

    // Need to check if the newton solver i converged.
    // If not, we half the stepsize and try again
    if (!step_ok){
      dt/=2.0;
      justrefined=true;
#ifdef DEBUG
      printf("First node failed, new dt = %1.6e\n",dt);
      newtonIter2.push_back(0);
      newtonAccepted1.push_back(0);
      newtonAccepted2.push_back(0);
#endif
      continue;
    }  
#ifdef DEBUG
    else{
      newtonAccepted1.push_back(1);
    }
#endif

    for (i = 0; i <ode->size(); ++i)
      prev[i]=y[i]+z1[i];

    ode->eval(prev,t+c1*dt,f1);

    for (i = 0; i <ode->size(); ++i)
      prev[i] = a21*f1[i];

    //Use z1 as initial guess for z2:
    for (i = 0; i <ode->size(); ++i)
      z2[i] = z1[i];

    step_ok = NewtonSolve(z2,prev,y,t+c2*dt,dt,a22);    
#ifdef DEBUG
    newtonIter2.push_back(newtonits);
#endif


    if (!step_ok){
      dt/=2.0;
      justrefined=true;
#ifdef DEBUG
      printf("Second node failed, new dt = %1.6e\n",dt);
      newtonAccepted2.push_back(0);
#endif
      continue;
    }else{
      t+=dt;
#ifdef DEBUG
      newtonAccepted2.push_back(1);
#endif
      if (fabs(t-t_end)<eps)
        done = true;
      else{
        // If the solver has refined, we do not allowe it to double its timestep for anoter step
        if (justrefined){
          justrefined=false;
        }else{
          double tmp = 2.0*dt;
          if (fabs(ldt-tmp)<eps)
            dt=ldt;
          else
            dt=tmp;
        }
        if ((t+ldt)>t_end)
          dt=t_end-t;
      }
    }


    for (i = 0; i < ode->size(); ++i)
      y[i] += (d1*z1[i]+d2*z2[i]);

    //printf("End of time loop\n");

  }  
#ifdef DEBUG
  printf("SDIRK2O3 done with comp_jac = %d and rejected = %d at t=%1.2e in %ld steps\n",jac_comp,rejects,t,num_tsteps);
#endif
}
