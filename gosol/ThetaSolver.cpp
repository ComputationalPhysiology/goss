#include "ThetaSolver.h"
//#include <iostream.h>
#include <stdio.h>

using namespace gosol;

//-----------------------------------------------------------------------------
ThetaSolver:: ThetaSolver (ODE* ode_, double _ldt)
  //-----------------------------------------------------------------------------
{ 
  attach(ode_);
  ldt = _ldt; 
} 

//-----------------------------------------------------------------------------
ThetaSolver:: ThetaSolver(double _ldt)
  //-----------------------------------------------------------------------------
{
  ldt = _ldt;
}

//-----------------------------------------------------------------------------
ThetaSolver:: ~ThetaSolver ()
  //-----------------------------------------------------------------------------
{
  free(z1); 
  free(ft1);
}


//-----------------------------------------------------------------------------
void  ThetaSolver:: attach(ODE* ode_)
  //-----------------------------------------------------------------------------
{
  ode = ode_;
  //  ode->attach(this);
  justrefined = false;
  num_tsteps = 0;
  theta = 0.5;
  init();

  stages = 1;
  z1 =  static_cast<double*>(malloc(sizeof*z1*ode->size())); 
  ft1 = static_cast<double*>(malloc(sizeof*ft1*ode->size())); 
}


//-----------------------------------------------------------------------------
void ThetaSolver:: forward(double* y, double t, double dt)
  //-----------------------------------------------------------------------------
{
  int i;
  double t_end = t+dt;
  dt = ldt;
  for (i = 0; i < ode->size(); ++i)
    prev[i] = 0.0;

  bool step_ok, done=false;
  double eps = 1e-14;// A way to chech if we are at t_end.

  while (!done){
    num_tsteps+=1;
    if (recompute_jacobian){
      computeJacobian(t,y);
      mult(-dt*(1-theta),jac);
      addIdentity(jac);
      factLU(jac);
      jac_comp+=1;
    }

    //Use 0.0 z1:
    for (i = 0; i <ode->size(); ++i)
      z1[i] = 0.0;

    ode->eval(y,t,ft1);
    for (i = 0; i <ode->size(); ++i){
      prev[i] = theta*ft1[i];
    }
    step_ok = NewtonSolve(z1,prev,y,t+dt,dt,theta);    
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
        y[i] += z1[i];
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
  printf("ThetaSolver done with comp_jac = %d and rejected = %d at t=%1.2e in %ld steps\n",jac_comp,rejects,t, num_tsteps);
#endif

}



