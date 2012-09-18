#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include "types.h"
#include "RKF32.h"

using namespace goss;

//-----------------------------------------------------------------------------
RKF32::RKF32() 
  : AdaptiveExplicitSolver(), 
    nfevals(0), ndtsa(0), ndtsr(0),
    a21(1.0/2.0), 
    a32(3.0/4.0), 
    b1(2.0/9.0), 
    b2(1.0/3.0), 
    b3(4.0/9.0), 
    bh1(7.0/24.0), 
    bh2(1.0/4.0), 
    bh3(1.0/3.0), 
    bh4(1.0/8.0), 
    d1(b1-bh1), 
    d2(b2-bh2), 
    d3(b3-bh3), 
    d4(-bh4), 
    c2(1.0/2.0), 
    c3(3.0/4.0), 
    nbytes(0),
    ki(0), k1(0), k2(0), k3(0), k4(0), yn(0), e(0)
{
  _iord = 3;
}
//-----------------------------------------------------------------------------
RKF32::RKF32(double ldt) 
  : AdaptiveExplicitSolver(ldt), 
    nfevals(0), ndtsa(0), ndtsr(0),
    a21(1.0/2.0), 
    a32(3.0/4.0), 
    b1(2.0/9.0), 
    b2(1.0/3.0), 
    b3(4.0/9.0), 
    bh1(7.0/24.0), 
    bh2(1.0/4.0), 
    bh3(1.0/3.0), 
    bh4(1.0/8.0), 
    d1(b1-bh1), 
    d2(b2-bh2), 
    d3(b3-bh3), 
    d4(-bh4), 
    c2(1.0/2.0), 
    c3(3.0/4.0), 
    nbytes(0),
    ki(0), k1(0), k2(0), k3(0), k4(0), yn(0), e(0)
{
  _iord = 3;
}
//-----------------------------------------------------------------------------
RKF32::RKF32(ODE* ode, double ldt) 
  : AdaptiveExplicitSolver(ldt, 0.0), 
    nfevals(0), ndtsa(0), ndtsr(0),
    a21(1.0/2.0), 
    a32(3.0/4.0), 
    b1(2.0/9.0), 
    b2(1.0/3.0), 
    b3(4.0/9.0), 
    bh1(7.0/24.0), 
    bh2(1.0/4.0), 
    bh3(1.0/3.0), 
    bh4(1.0/8.0), 
    d1(b1-bh1), 
    d2(b2-bh2), 
    d3(b3-bh3), 
    d4(-bh4), 
    c2(1.0/2.0), 
    c3(3.0/4.0), 
    nbytes(0),
    ki(0), k1(0), k2(0), k3(0), k4(0), yn(0), e(0)
{ 
  _iord = 3;
  attach(ode);
}
//-----------------------------------------------------------------------------
RKF32::~RKF32() 
{
  if (ki) delete[] ki;
  if (k1) delete[] k1;
  if (k2) delete[] k2;
  if (k3) delete[] k3;
  if (k4) delete[] k4;
  if (yn) delete[] yn;
  if (e)  delete[] e;
}

//-----------------------------------------------------------------------------
void RKF32::attach(ODE* ode)
{
  // Attach ode using base class. 
  // NOTE: This will trigger call to reset
  ODESolver::attach(ode);

  // Delete memory if excist
  if (ki) delete[] ki;
  if (k1) delete[] k1;
  if (k2) delete[] k2;
  if (k3) delete[] k3;
  if (k4) delete[] k4;
  if (yn) delete[] yn;
  if (e)  delete[] e;

  // Initilize RK increments
  ki = new double[ode_size()]; 
  k1 = new double[ode_size()];
  k2 = new double[ode_size()];
  k3 = new double[ode_size()];
  k4 = new double[ode_size()];
  yn = new double[ode_size()];
  e  = new double[ode_size()];

  nbytes  = ode_size()*sizeof(double);
}
//-----------------------------------------------------------------------------
void RKF32::reset()
{
  first = true;

  nfevals = 0;
  ndtsa   = 0;
  ndtsr   = 0;
  num_accepted = 0;
  num_rejected = 0;

  // Reset base classes
  AdaptiveExplicitSolver::reset();
}
//-----------------------------------------------------------------------------
void RKF32::forward(double* y, double t, double interval) 
{
  uint i;

  // We swap the result vektor if a timestep is accepted. We therefore need 
  // to store the pointer to the initial y-vector in order to ensure that 
  // this memory segment contains the final result when the end is reached ...
  double* ret_ptr = y;
  double* swap;

  // End of interval
  const double t_end = t + interval;

  reached_tend = false;

  // Local time It is only used in 
  _t = t;

  // FIXME: first i always true
  if (first) 
  {
    _ode->eval(y, t, k1);
    nfevals += 1;
  }

  // Set initial time step
  if (_ldt < 0) 
  {
    _dt = dtinit(t, y, yn, k1, k2, _iord);
    nfevals += 1;
  } 
  else 
  {
    _dt = _ldt;
  }

#ifdef DEBUG
  // log data
  dt_v.push_back(_dt);
#endif

  while (!reached_tend)
  {
    
    for (i = 0; i < ode_size(); ++i)
      ki[i] = y[i] + _dt*a21*k1[i];

    _ode->eval(ki, t + c2*_dt, k2);
    
    for (i = 0; i < ode_size(); ++i)
      ki[i] = y[i] + _dt*a32*k2[i];

    _ode->eval(ki, t + c3*_dt, k3);

    // We assemble the new y
    for (i = 0; i < ode_size(); ++i)
      yn[i] = y[i] + _dt*(b1*k1[i] + b2*k2[i] + b3*k3[i]);
    //yn[i] = y[i] + _dt*(bh1*k1[i]+bh2*k2[i]+bh3*k3[i]+bh4*k4[i]);

    // We compute the first quadrature node for the next iteration (FSAL)
    _ode->eval(yn, t + _dt, k4);
    nfevals += 3;

    // We compute the error vector
    for (i=0; i < ode_size(); ++i)
      e[i] = _dt*(d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i]);

    // Compute new time step and check if it is rejected
    new_time_step(y, yn, e, t_end);

#ifdef DEBUG
    log_data(_dt, step_accepted);
#endif

    // If step
    if (step_accepted)
    {
      ndtsa += 1;
      swap = y;
      y    = yn;
      yn   = swap;
      swap = k1;
      k1   = k4;
      k4   = swap;
#ifdef DEBUG
      if (single_step_mode)
      {
        if (ret_ptr != y)
	{
          memcpy(ret_ptr, y, nbytes);
          yn = y;
        }
        swap = 0;
        return;
      }
#endif
    } 
    else
      ndtsr += 1;

  }

  // This is a copy to ensure that the input Ptr contains the final solution
  // This can probably be done in a more elegant way
  if (ret_ptr != y)
  {
    memcpy(ret_ptr, y, nbytes);
    yn = y;
  }

#ifdef DEBUG
  dt_v.pop_back();
#endif
}
//-----------------------------------------------------------------------------
#ifdef DEBUG
void RKF32::log_data(double dt, bool accepted)
{
  dt_v.push_back(dt);
  accept_v.push_back(accepted);
}
//-----------------------------------------------------------------------------
void RKF32::dt_vector(DoubleVector *res)
{
  res->n    = dt_v.size();
  res->data = new double[dt_v.size()]; 
  for(uint i = 0; i < dt_v.size(); ++i)
    res->data[i] = dt_v[i];
}

//-----------------------------------------------------------------------------
void RKF32::accepted_vector(DoubleVector *res)
{
  res->n    = accept_v.size();
  res->data = new double[accept_v.size()]; 
  for(uint i = 0; i < accept_v.size(); ++i)
    res->data[i] = float(accept_v[i]);
}

#else
//-----------------------------------------------------------------------------
void RKF32::log_data(double, bool)
{
  std::cout << "DEBUG OFF!" << std::endl; 
}

//-----------------------------------------------------------------------------
void RKF32::dt_vector(DoubleVector*)
{
  std::cout << "DEBUG OFF!" << std::endl; 
}
//-----------------------------------------------------------------------------
void RKF32::accepted_vector(DoubleVector*)
{
  std::cout << "DEBUG OFF!" << std::endl; 
}
//-----------------------------------------------------------------------------
#endif
