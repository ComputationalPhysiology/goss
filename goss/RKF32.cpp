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
RKF32::RKF32(const RKF32& solver) 
  : AdaptiveExplicitSolver(solver), 
    nfevals(solver.nfevals), ndtsa(solver.ndtsa), ndtsr(solver.ndtsr),
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
    nbytes(solver.num_states()*sizeof(double)),
    ki(new double[solver.num_states()]), k1(new double[solver.num_states()]), 
    k2(new double[solver.num_states()]), k3(new double[solver.num_states()]), 
    k4(new double[solver.num_states()]), yn(new double[solver.num_states()]), 
    e(new double[solver.num_states()])
{ 
  // Do nothing
}
//-----------------------------------------------------------------------------
RKF32::~RKF32() 
{
  // Do nothing
}

//-----------------------------------------------------------------------------
void RKF32::attach(ODE* ode)
{
  // Attach ode using base class. 
  // NOTE: This will trigger call to reset
  ODESolver::attach(ode);

  // Initilize RK increments
  ki.reset(new double[num_states()]); 
  k1.reset(new double[num_states()]);
  k2.reset(new double[num_states()]);
  k3.reset(new double[num_states()]);
  k4.reset(new double[num_states()]);
  yn.reset(new double[num_states()]);
  e .reset(new double[num_states()]);

  nbytes  = num_states()*sizeof(double);
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

  // We swap the result vector if a timestep is accepted. We therefore need 
  // to store the pointer to the initial y-vector in order to ensure that 
  // this memory segment contains the final result when the end is reached ...
  double* ret_ptr = y;
  double* swap, *yn0;

  // Use the raw pointer instead of the scoped_array, enabling pointer swap
  yn0 = yn.get();

  // End of interval
  const double t_end = t + interval;

  reached_tend = false;

  // Local time It is only used in 
  _t = t;

  // FIXME: first i always true
  if (first) 
  {
    _ode->eval(y, t, k1.get());
    nfevals += 1;
  }

  // Set initial time step
  if (_ldt < 0) 
  {
    _dt = dtinit(t, y, yn0, k1.get(), k2.get(), _iord);
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
    
    for (i = 0; i < num_states(); ++i)
      ki[i] = y[i] + _dt*a21*k1[i];

    _ode->eval(ki.get(), t + c2*_dt, k2.get());
    
    for (i = 0; i < num_states(); ++i)
      ki[i] = y[i] + _dt*a32*k2[i];

    _ode->eval(ki.get(), t + c3*_dt, k3.get());

    // We assemble the new y
    for (i = 0; i < num_states(); ++i)
      yn0[i] = y[i] + _dt*(b1*k1[i] + b2*k2[i] + b3*k3[i]);
    //yn[i] = y[i] + _dt*(bh1*k1[i]+bh2*k2[i]+bh3*k3[i]+bh4*k4[i]);

    // We compute the first quadrature node for the next iteration (FSAL)
    _ode->eval(yn0, t + _dt, k4.get());
    nfevals += 3;

    // We compute the error vector
    for (i=0; i < num_states(); ++i)
      e[i] = _dt*(d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i]);

    // Compute new time step and check if it is rejected
    new_time_step(y, yn0, e.get(), t_end);

#ifdef DEBUG
    log_data(_dt, step_accepted);
#endif

    // If step
    if (step_accepted)
    {
      ndtsa += 1;
      swap = y;
      y = yn0;
      yn0 = swap;
      k4.swap(k1);
#ifdef DEBUG
      if (single_step_mode)
      {
        if (ret_ptr != y)
	{
          memcpy(ret_ptr, y, nbytes);
          yn0 = y;
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
    yn0 = y;
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
  res->data.reset(new double[dt_v.size()]);
  for(uint i = 0; i < dt_v.size(); ++i)
    res->data[i] = dt_v[i];
}

//-----------------------------------------------------------------------------
void RKF32::accepted_vector(DoubleVector *res)
{
  res->n    = accept_v.size();
  res->data.reset(new double[accept_v.size()]); 
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
