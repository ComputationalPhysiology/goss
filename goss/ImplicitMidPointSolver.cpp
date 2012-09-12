#include "ImplicitMidPointSolver.h"
//#include <iostream.h>
#include <stdio.h>

using namespace goss;

void  ImplicitMidPointSolver:: attach(ODE* ode_){
  ode = ode_;
  justrefined=false;
  num_tsteps=0;
  init();
  stages = 1;
  z1 = static_cast<double*>(malloc(sizeof(double)*ode->size())); 
}


void ImplicitMidPointSolver:: forward(double* y, double t, double dt) {
  int i;
  double t_end = t + dt;
  dt = ldt;
  for (i = 0; i < ode->size(); ++i)
    prev[i] = 0.0;

  bool step_ok, done=false;
  double eps = 1e-14;// A way to chech if we are at t_end.
  //printf("t_start=%1.2e, end t_end=%1.2e, and dt=%1.2e\n",t,t_end,dt);

  while (!done){
    num_tsteps+=1;
    if (recompute_jacobian){
      computeJacobian(t,y);
      mult(-dt*0.5,jac);
      addIdentity(jac);
      factLU(jac);
      jac_comp+=1;
    }

    //Use 0.0 z1:
    for (i = 0; i <ode->size(); ++i)
      z1[i] = 0.0;

    step_ok = NewtonSolve(z1,prev,y,t+0.5*dt,dt,0.5);    
#ifdef DEBUG
    newtonIter1.push_back(newtonits);
    dt_v.push_back(dt);
#endif
    if (step_ok){
      t+=dt;
      if (fabs(t-t_end)<eps){
        done = true;
        //printf("Done at t=%1.4e\n",t);
      }else{
        // If the solver has refined, we do not allowe it to double its timestep for anoter step
        if (justrefined){
          justrefined=false;
        }else{
          double tmp = 2.0*dt;
          //if (fabs(ldt-tmp)<eps)
          if (tmp>=ldt)
            dt=ldt;
          else
            dt=tmp;
        }
        if ((t+dt)>t_end)
          dt=t_end-t;
      }
      for (i = 0; i < ode->size(); ++i)
        y[i] += 2*z1[i]; //Remember that in this formulation, 
      // y_n=y_n-1+sum d_i z_i, where d=bA^-1
#ifdef DEBUG
      newtonAccepted1.push_back(1);
#endif    
    }else{
      dt/=2.0;
      justrefined=true;
#ifdef DEBUG
      newtonAccepted1.push_back(0);
#endif    
    }
  }
#ifdef DEBUG
  printf("ImplicitMidPointSolver done with comp_jac = %d and rejected = %d at t=%1.2e in %ld steps\n",jac_comp,rejects,t, num_tsteps);
#endif

}



