//
// $Id$
// 

#include "ESDIRK4O32.h"
#include <stdio.h>

using namespace gosol;

void  ESDIRK4O32:: attach(ODE* ode_){

  ode = ode_;
  //ode->attach(this);
  int system_size = ode->size();
  //init();

  gamma = 0.43586652150845899942;
  //newton_tol = 1e-5;

  a21 = gamma;
  a22 = gamma;
  a31 = (-4*gamma*gamma+6*gamma-1)/(4*gamma);
  a32 = (-2*gamma+1)/(4*gamma);
  a33 = gamma;
  a41 = (6*gamma-1)/(12*gamma);
  a42 = -1/((24*gamma-12)*gamma);
  a43 = (-6*gamma*gamma+6*gamma-1)/(6*gamma-3);
  a44 = gamma;
  bh1 = a31;
  bh2 = a32;
  bh3 = a33;
  b1  = a41;
  b2  = a42;
  b3  = a43;
  b4  = a44;
  c2  = 2.0*gamma;
  c3  = 1.0;
  c4  = 1.0;

  z1 = static_cast<double*>(malloc(sizeof(double)*system_size)); 
  z2 = static_cast<double*>(malloc(sizeof(double)*system_size)); 
  z3 = static_cast<double*>(malloc(sizeof(double)*system_size)); 
  z4 = static_cast<double*>(malloc(sizeof(double)*system_size)); 
  yn = static_cast<double*>(malloc(sizeof(double)*system_size)); 
  yh = static_cast<double*>(malloc(sizeof(double)*system_size)); 
  nfevals = 0;
  ndtsa   = 0;
  ndtsr   = 0;
}


void ESDIRK4O32:: forward(double* y, double t_, double DT) {
  //NB sjekk definisjonen av prev vs hvilken dt som sendes til NewtonSolve!

  t=t_;
  int i;
  double t_end = t+DT;
  //if (ldt>0) {
  //    // Calculates a local time step that is <= ldt,
  //    // and that divides dt into equally sized steps:
  //
  //    n_steps = (int) ceil(dt/ldt -1E-12); 
  //    dt /= n_steps; //could have used new variable, but dt was available...
  //}


  dt=ldt;    
  bool step_ok; // done = false;
  reached_tend=false;

  retPtr = y;
  if (dt<0){
    ode->eval(y,t,z1);
    dt = dtinit(t, y, yn, z1, z2, iord);
    nfevals+=2;
  }else{
    dt=ldt;
  }


#ifdef DEBUG
  // log data
  dt_v.push_back(dt);
#endif

  while (!reached_tend){
    //Jacobian recomputed once for each local step
    if (recompute_jacobian){
      computeJacobian(t,y);
      nfevals += ode->size();
      //computeJacobian(t+c2*dt,y);
      mult(-dt*a22,jac);
      addIdentity(jac);
      factLU(jac);
      jac_comp+=1;
      recompute_jacobian = false;
      //printf("Recomputed jac at t=%1.4e with dt=%1.4e\n",t,dt);
    }

    // Computes the first explicit node
    ode->eval(y,t,z1);
    nfevals += 1;

    //printf("z2=");
    // Computes the second node, implicitly
    for (i = 0; i < ode->size(); ++i){
      prev[i] = a21*z1[i];
      z2[i]   = 0.0; //z1[i]-y[i]; // Ola comment: What happens here??
      //printf("%1.2e, ",z2[i]);
    }

    step_ok = NewtonSolve(z2,prev,y,t+c2*dt,dt,a22);    


    // Need to check if the newton solver is converged.
    // If not, we half the stepsize and try again
    if (!step_ok){
      recompute_jacobian=true;
      dt/=2.0;
#ifdef DEBUG
      logData(dt,false);
      printf("Second node failed, new dt = %1.6e\n",dt);
#endif
      continue;
    }

    //printf("z2=");
    for (i = 0; i < ode->size(); ++i){
      z3[i]=z2[i];
      z2[i]+=y[i];
      //printf("%1.2e, ",z2[i]);
    }
    //printf(" before eval\n");
    ode->eval(z2,t+c2*dt,f1);
    nfevals+=1;

    // Computes the third node, implicitly
    // Use pointer swap instead of copy!!
    //printf("z2=");
    for (i = 0; i <ode->size(); ++i){
      //printf("%1.2e, ",f1[i]);
      z2[i]=f1[i];
      prev[i] = a31*z1[i]+a32*z2[i];
    }
    //printf(" after eval\n");

    step_ok = NewtonSolve(z3,prev,y,t+c3*dt,dt,a33);

    // Need to check if the newton solver is converged.
    // If not, we half the stepsize and try again
    if (!step_ok){
      recompute_jacobian=true;
      dt/=2.0;
#ifdef DEBUG
      logData(dt,false);
      printf("Third node failed, new dt = %1.6e\n",dt);
#endif
      continue;
    }

    for (i = 0; i < ode->size(); ++i)
      z3[i]+=y[i];

    ode->eval(z3,t+c3*dt,f1);
    nfevals+=1;
    // Computes the error estimate
    // Use pointer swap instead of copy!!
    for (i = 0; i <ode->size(); ++i){
      z3[i]=f1[i];
      yh[i]=y[i]+dt*(bh1*z1[i]+bh2*z2[i]+bh3*z3[i]);
    }

    // Computes the fourth node, implicitly
    for (i = 0; i < ode->size(); ++i){
      z4[i]=yh[i]-y[i];
      prev[i] = a41*z1[i]+a42*z2[i]+a43*z3[i];
    }

    step_ok = NewtonSolve(z4,prev,y,t+c4*dt,dt,a44);

    // Need to check if the newton solver is converged.
    // If not, we half the stepsize and try again
    if (!step_ok){
      recompute_jacobian=true;
      dt/=2.0;
#ifdef DEBUG
      logData(dt,false);
      printf("Fourth node failed, new dt = %1.6e\n",dt);
#endif
      continue;
    } else {
      for (i = 0; i < ode->size(); ++i)
        z4[i]+=y[i];
      ode->eval(z4,t+c4*dt,f1);
      nfevals+=1;
      for (i = 0; i <ode->size(); ++i){
        //printf("i=%d, dt=%1.4e, b1=%1.4e, b2=%1.4e, b3=%1.4e, b4=%1.4e\n",i,dt,b1,b2,b3,b4);
        //printf("yn=%1.4e, y=%1.4e, z1=%1.4e, z2=%1.4e, z3=%1.4e, z4=%1.4e\n",yn[i],y[i],z1[i] ,z2[i] ,z3[i] ,f1[i]); 
        yn[i] = y[i]+dt*(b1*z1[i]+b2*z2[i]+b3*z3[i]+b4*f1[i]);
        yh[i]-=yn[i];//use this as the error vector
      }

      newTimeStep( y, yn, yh, t_end);//,blog);
      //recompute_jacobian=true;

      //done=blog[1];
#ifdef DEBUG
      logData(dt,step_accepted);
#endif
      if (step_accepted){
        swap = y;
        y    = yn;
        yn   = swap;
#ifdef DEBUG
        if (single_step_mode){
          if (retPtr!=y){
            for (i=0; i < ode->size(); ++i)
              retPtr[i] = y[i];
            yn = y;
          }
          swap=0;
          return;
        }
#endif
      }


    }



  }
  // This is a copy to ensure that the input Ptr contains the final solution
  // This can probably be done in a more elegant way
  if (retPtr!=y){
    for (i=0; i < ode->size(); ++i)
      retPtr[i] = y[i];
    yn = y;
    //printf("Accepted at t=%1.4e, with dt=%1.4e\n",t,dt);
  }

#ifdef DEBUG
  dt_v.pop_back();
  //cout << "accepted,rejected="<<log[0]<<","<<log[1]<<endl;
  printf("ESDIRK4O32 done with comp_jac = %d and rejected = %d\n",jac_comp,rejects);
#endif


  //printf("Jumped out of the time loop at t = %1.6e\n",t);
}





//
// $Log$
// 
