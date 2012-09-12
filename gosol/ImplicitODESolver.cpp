// This class needs some serious reconsidderation regarding memory use!!


#include "ImplicitODESolver.h"
//#include <iostream.h>
#include <stdio.h>


using namespace gosol;

ImplicitODESolver:: ImplicitODESolver()
  : ODESolver(0.0, 0.0), jac(0), f1(0), f2(0), yz(0), b(0), dz(0), prev(0),
    newton_tol(0.0), eta(1.0), kappa(0.1), stages(0), newtonits(0), maxits(0), rejects(0),
    jac_comp(0), min_dt(0.0), recompute_jacobian(true)
{
  // Do nothing
}


ImplicitODESolver:: ImplicitODESolver(double ldt)
  : ODESolver(ldt), jac(0), f1(0), f2(0), yz(0), b(0), dz(0), prev(0),
    newton_tol(0.0), eta(1.0), kappa(0.1), stages(0), newtonits(0), maxits(0), rejects(0),
    jac_comp(0), min_dt(0.0), recompute_jacobian(true)
{
  // Do nothing
}

void ImplicitODESolver:: init(){
  maxits = 10;
  rejects = 0;
  jac_comp = 0;
  eta      = 1e-10;
  kappa = 1e-1;
  recompute_jacobian = true;//We need to compute the Jacobian the first time



  b    = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  dz   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  prev = static_cast<double*>(malloc(sizeof(double)*ode->size()));  

  yz   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  f1   = static_cast<double*>(malloc(sizeof(double)*ode->size()));
  f2   = static_cast<double*>(malloc(sizeof(double)*ode->size()));

  jac  = static_cast<double**>(malloc(sizeof(double*)*ode->size()));
  jac_size = ode->size();
  for (int i = 0; i < ode->size(); ++i)
    jac[i] = static_cast<double*>(malloc(sizeof(double)*ode->size()));

  setNewtonTol();

}

void ImplicitODESolver:: computeJacobian(double t, double* y){

  int i,j;
  double max,ysafe,delta;
  ode->eval(y,t, f1);
  for (i=0;i<ode->size();++i){
    ysafe = y[i];
    max=((1e-5>fabs(ysafe))?1e-5:fabs(ysafe));
    delta=sqrt(1e-15*max);
    y[i] += delta;
    ode->eval(y, t, f2);
    for (j=0;j<ode->size();++j){
      jac[j][i]=(f2[j]-f1[j])/delta;
    }
    y[i]=ysafe;
  } 
}

void ImplicitODESolver:: mult(double fact, double** matrix){
  int i,j;
  for (i = 0; i < ode->size(); ++i)
    for (j = 0; j < ode->size(); ++j)
      matrix[i][j] *= fact;
}

void ImplicitODESolver:: addIdentity(double** matrix){
  int i;
  for (i = 0; i < ode->size(); ++i)
    matrix[i][i] += 1;
}


void ImplicitODESolver:: factLU(double** mat){
  int i, k, r;
  double sum;


  for (k = 1; k < ode->size(); k++){
    for (i = 0; i <= k-1; ++i){
      sum = 0.0;
      for (r = 0; r <= i-1; r++)
        sum += mat[i][r]*mat[r][k];
      mat[i][k] -=sum;
      sum = 0.0;
      for (r = 0; r <= i-1; r++)
        sum += mat[k][r]*mat[r][i];
      mat[k][i] = (mat[k][i]-sum)/mat[i][i];
    }
    sum = 0;
    for (r = 0; r <= k-1; r++)
      sum += mat[k][r]*mat[r][k];
    mat[k][k] -= sum;
  }


}

void ImplicitODESolver::forwBackLU(const double * const * mat, double* b, double* x){
  //solves Ax = b with forward backward substitution, provided that 
  //A is already L1U factorized

  double sum;

  int i,j;

  x[0] = b[0];

  for (i = 1; i < ode->size(); ++i){
    sum = 0;
    for (j = 0; j <= i-1; ++j)
      sum = sum + mat[i][j]*x[j];
    x[i] = b[i] -sum;
  }
  x[ode->size()-1] = x[ode->size()-1]/mat[ode->size()-1][ode->size()-1];


  for (i = ode->size()-2; i >=0; i--){
    sum = 0;
    for (j = i+1; j < ode->size(); ++j)
      sum = sum +mat[i][j]*x[j];
    x[i] = (x[i]-sum)/mat[i][i];
  }

}


bool ImplicitODESolver::NewtonSolve(double* z, double* prev, double* y0, double t, double dt, double alpha){
  int i;
  bool step_ok = true, converged = false;
  newtonits = 0;
  double Ntheta, z_norm, prev_norm;
  recompute_jacobian=false;

  do{
    for (i = 0; i < ode->size(); ++i)
      yz[i]=y0[i]+z[i];
    ode->eval(yz,t,f1);
    for (i = 0; i < ode->size(); ++i)
      b[i] = -z[i]+dt*(prev[i]+alpha*f1[i]);

    forwBackLU(jac,b,dz);
    z_norm = norm(dz);



    if (newtonits > 0) {
      Ntheta = z_norm/prev_norm;
      if (Ntheta<1e-3)
        recompute_jacobian=false;
      else
        recompute_jacobian=true;
      if (Ntheta >1){
#ifdef DEBUG
        printf("Newton solver diverges with Ntheta = %1.2e, reduces time step\n", Ntheta);
#endif
        rejects ++;
        step_ok = false;
        //return step_ok;
        break;
      }
      if (z_norm>(kappa*newton_tol*(1-Ntheta)/pow(Ntheta,maxits-newtonits))){
#ifdef DEBUG
        printf("Newton solver converges to slow with theta = %1.2e at iteration %d, reduces time step\n", Ntheta,newtonits);
#endif
        rejects ++;
        step_ok = false;
        recompute_jacobian=true;
        break;
      }
      eta = Ntheta/(1.0-Ntheta);
    }else{
      eta = (eta > 1e-15 ? eta : 1e-15);
      eta = pow(eta,0.8);
    }

    if (newtonits > maxits && !converged){
#ifdef DEBUG
      printf("Not converged in %d iterations. Reduces time step.\n",maxits);
#endif
      rejects ++;
      step_ok = false;
      //return step_ok;
      break;
    }
    for (i = 0; i <ode->size(); ++i)
      z[i] += dz[i];

    prev_norm = z_norm;
    newtonits++;
  } while (eta*z_norm <= kappa*newton_tol); 


  return step_ok;
}





double ImplicitODESolver:: norm(double* vec)
{
  double l2_norm;

  for (int i = 0; i < ode->size(); ++i)
    l2_norm += vec[i]*vec[i];

  l2_norm = sqrt(l2_norm);
  return l2_norm;
}
