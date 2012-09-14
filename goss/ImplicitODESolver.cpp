// FIXME: This class needs some serious reconsideration regarding memory use!!

#include "ImplicitODESolver.h"
//#include <iostream.h>
#include <stdio.h>
#include <cmath>

using namespace goss;

//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver()
  : ODESolver(0.0, 0.0), jac(0), f1(0), f2(0), yz(0), _b(0), dz(0), _prev(0),
    _newton_tol(1.e-5), eta(1e-10), kappa(0.1), stages(0), newtonits(0), 
    maxits(10), rejects(0), jac_comp(0), min_dt(0.0), recompute_jacobian(true)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitODESolver::ImplicitODESolver(double ldt)
  : ODESolver(ldt), jac(0), f1(0), f2(0), yz(0), _b(0), dz(0), _prev(0),
    _newton_tol(1.e-5), eta(1e-10), kappa(0.1), stages(0), newtonits(0), 
    maxits(10), rejects(0), jac_comp(0), min_dt(0.0), recompute_jacobian(true)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ImplicitODESolver::~ImplicitODESolver ()
{
  // Destroy jacobian
  if (jac) 
  {
    for (uint i=0; i < jac_size; ++i)
      delete[] jac[i];
  }
  
  if (jac)  delete[] jac;
  if (f1)   delete[] f1;
  if (f2)   delete[] f2;
  if (yz)   delete[] yz;
  if (_b)    delete[] _b; 
  if (dz)   delete[] dz;
  if (_prev) delete[] _prev;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::init()
{
  rejects = 0;
  jac_comp = 0;

  // We need to compute the Jacobian the first time
  recompute_jacobian = true;

  _b   = new double[_ode->size()]; 
  dz   = new double[_ode->size()];
  _prev = new double[_ode->size()];

  yz   = new double[_ode->size()];
  f1   = new double[_ode->size()];
  f2   = new double[_ode->size()];

  // Init jacobian
  jac_size = _ode->size();
  jac  = new double*[_ode->size()];
  for (uint i = 0; i < _ode->size(); ++i)
    jac[i] = new double[_ode->size()];

  // Newton tolerance
  _newton_tol = 1.e-5;
  
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::compute_jacobian(double t, double* y)
{
  uint i, j;
  double max, ysafe, delta;
  _ode->eval(y, t, f1);
  
  for (i = 0; i < _ode->size(); ++i)
  {
    ysafe = y[i];
    max = 1e-5 > std::fabs(ysafe) ? 1e-5 : std::fabs(ysafe);
    delta = std::sqrt(1e-15*max);
    y[i] += delta;
    _ode->eval(y, t, f2);
    
    for (j=0;j<_ode->size();++j)
      jac[j][i]=(f2[j]-f1[j])/delta;
    
    y[i]=ysafe;
  } 
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::mult(double fact, double** matrix)
{
  for (uint i = 0; i < _ode->size(); ++i)
    for (uint j = 0; j < _ode->size(); ++j)
      matrix[i][j] *= fact;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::add_identity(double** matrix)
{
  for (uint i = 0; i < _ode->size(); ++i)
    matrix[i][i] += 1;
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::lu_factorize(double** mat)
{
  double sum;

  for (uint k = 1; k < _ode->size(); k++)
  {

    for (uint i = 0; i <= k-1; ++i)
    {
      sum = 0.0;
      for (uint r = 0; r <= i-1; r++)
        sum += mat[i][r]*mat[r][k];
    
      mat[i][k] -=sum;
      sum = 0.0;

      for (uint r = 0; r <= i-1; r++)
        sum += mat[k][r]*mat[r][i];
    
      mat[k][i] = (mat[k][i]-sum)/mat[i][i];

    }

    sum = 0.0;
    for (uint r = 0; r <= k-1; r++)
      sum += mat[k][r]*mat[r][k];

    mat[k][k] -= sum;

  }
}
//-----------------------------------------------------------------------------
void ImplicitODESolver::forward_backward_subst(const double * const * mat, 
					       double* b, double* x)
{
  // solves Ax = b with forward backward substitution, provided that 
  // A is already L1U factorized

  double sum;

  x[0] = b[0];

  for (uint i = 1; i < _ode->size(); ++i)
  {
    sum = 0.0;
    for (uint j = 0; j <= i-1; ++j)
      sum = sum + mat[i][j]*x[j];

    x[i] = b[i] -sum;
  }

  x[_ode->size()-1] = x[_ode->size()-1]/mat[_ode->size()-1][_ode->size()-1];

  for (int i = _ode->size()-2; i >= 0; i--)
  {
    sum = 0;
    for (uint j = i + 1; j < _ode->size(); ++j)
      sum = sum +mat[i][j]*x[j];
  
    x[i] = (x[i]-sum)/mat[i][i];
  }
}
//-----------------------------------------------------------------------------
bool ImplicitODESolver::newton_solve(double* z, double* prev, double* y0, double t,
				     double dt, double alpha)
{
  uint i;
  bool step_ok = true, converged = false;
  newtonits = 0;
  double Ntheta, z_norm, prev_norm = 1.0;
  recompute_jacobian = false;

  do
  {
    for (i = 0; i < _ode->size(); ++i)
      yz[i] = y0[i] + z[i];

    _ode->eval(yz,t,f1);
    for (i = 0; i < _ode->size(); ++i)
      _b[i] = -z[i] + dt*(prev[i] + alpha*f1[i]);

    forward_backward_subst(jac, _b, dz);
    z_norm = norm(dz);

    if (newtonits > 0) 
    {
      Ntheta = z_norm/prev_norm;

      if (Ntheta<1e-3)
        recompute_jacobian = false;
      else
        recompute_jacobian = true;
    
      if (Ntheta >1)
      {
#ifdef DEBUG
        printf("Newton solver diverges with Ntheta = %1.2e, reduces time step\n", Ntheta);
#endif
        rejects ++;
        step_ok = false;
        //return step_ok;
        break;
      }
      
      if (z_norm > (kappa*_newton_tol*(1 - Ntheta)/std::pow(Ntheta, maxits - newtonits)))
      {
#ifdef DEBUG
        printf("Newton solver converges to slow with theta = %1.2e at "
	       "iteration %d, reduces time step\n", Ntheta, newtonits);
#endif
        rejects ++;
        step_ok = false;
        recompute_jacobian=true;
        break;
      }
      
      eta = Ntheta/(1.0-Ntheta);
    }
    
    // newtonits == 0
    else
    {
      eta = eta > 1e-15 ? eta : 1e-15;
      eta = std::pow(eta, 0.8);
    }

    // No convergence
    if (newtonits > maxits && !converged)
    {
#ifdef DEBUG
      printf("Not converged in %d iterations. Reduces time step.\n",maxits);
#endif
      rejects ++;
      step_ok = false;
      //return step_ok;
      break;
    }
    
    for (i = 0; i <_ode->size(); ++i)
      z[i] += dz[i];

    prev_norm = z_norm;
    newtonits++;

  } while (eta*z_norm <= kappa*_newton_tol); 

  return step_ok;
}
//-----------------------------------------------------------------------------
double ImplicitODESolver::norm(double* vec)
{
  double l2_norm = 0;

  for (uint i = 0; i < _ode->size(); ++i)
    l2_norm += vec[i]*vec[i];

  l2_norm = std::sqrt(l2_norm);
  return l2_norm;
}
//-----------------------------------------------------------------------------
