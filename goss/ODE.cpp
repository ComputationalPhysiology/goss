#include <iostream>
#include <stdexcept>
#include <cmath>

#include "Timer.h"
#include "log.h"
#include "ODE.h"

using namespace goss;

//-----------------------------------------------------------------------------
ODE::ODE(uint num_states_) : 
  _num_states(num_states_), _f1(num_states_), _f2(num_states_)
{ 
} 
//-----------------------------------------------------------------------------
double ODE::eval(uint idx, const double* states, double time)
{ 
  Timer timer_("Componentwise evaluation of rhs");

  warning("Warning: Calling base class ODE::eval component wise. "\
	  "This is very slow.");

  if (idx >= _num_states)
    error("Index out of range");

  eval(states, time, &_f1[0]);
  
  const double ret = _f1[idx];

  return ret;
}
//-----------------------------------------------------------------------------
void ODE::compute_jacobian(double* states, double time, double* jac)
{

  Timer timer_("Jacobian computation");

  uint i, j;
  double max, ysafe, delta;
  eval(states, time, &_f1[0]);
  
  for (i = 0; i < _num_states; ++i)
  {
    ysafe = states[i];
    max = 1e-5 > std::fabs(ysafe) ? 1e-5 : std::fabs(ysafe);
    delta = std::sqrt(1e-15*max);
    states[i] += delta;
    eval(states, time, &_f2[0]);
    
    for (j = 0; j < _num_states; ++j)
      jac[j*_num_states+i] = (_f2[j] - _f1[j])/delta;

    states[i] = ysafe;
  } 
}
//-----------------------------------------------------------------------------
void ODE::lu_factorize(double* mat) const
{

  Timer timer_("Factorizing jacobian");

  //std::cout << "Calling base class ODE::lu_factorize_subst." << std::endl;
  double sum;
  int i, k, r;
  int lnum_states = _num_states;

  for (k = 1; k < lnum_states; k++)
  {

    for (i = 0; i <= k-1; ++i)
    {
      sum = 0.0;
      for (r = 0; r <= i-1; r++)
        sum += mat[i*lnum_states+r]*mat[r*lnum_states+k];

      mat[i*lnum_states+k] -=sum;
      sum = 0.0;

      for (r = 0; r <= i-1; r++)
        sum += mat[k*lnum_states+r]*mat[r*lnum_states+i];
    
      mat[k*lnum_states+i] = (mat[k*lnum_states+i]-sum)/mat[i*lnum_states+i];

    }

    sum = 0.0;
    for (r = 0; r <= k-1; r++)
      sum += mat[k*lnum_states+r]*mat[r*lnum_states+k];

    mat[k*lnum_states+k] -= sum;

  }
}
//-----------------------------------------------------------------------------
void ODE::forward_backward_subst(const double* mat, const double* b, double* dx) const
{
  // solves Ax = b with forward backward substitution, provided that 
  // A is already LU factorized
  //std::cout << "Calling base class ODE::forward_backward_subst." << std::endl;

  Timer timer_("Forward backward substitution");

  double sum;

  dx[0] = b[0];

  for (uint i = 1; i < _num_states; ++i)
  {
    sum = 0.0;
    for (uint j = 0; j <= i-1; ++j)
      sum = sum + mat[i*_num_states+j]*dx[j];

    dx[i] = b[i] -sum;
  }

  const uint num_minus_1 = _num_states-1;
  dx[num_minus_1] = dx[num_minus_1]/mat[num_minus_1*_num_states+num_minus_1];

  for (int i = _num_states - 2; i >= 0; i--)
  {
    sum = 0;
    for (uint j = i + 1; j < _num_states; ++j)
      sum = sum +mat[i*_num_states+j]*dx[j];
  
    dx[i] = (dx[i]-sum)/mat[i*_num_states+i];
  }
}
//-----------------------------------------------------------------------------
void ODE::linear_terms(uint*) const
{
  error("ODE::linear_terms must be implement in a subclass");
}
//-----------------------------------------------------------------------------
void ODE::linear_derivatives(const double*, double, double*) const
{
  error("ODE::linear_derivatives must be implemented in a subclass");
}
//-----------------------------------------------------------------------------
