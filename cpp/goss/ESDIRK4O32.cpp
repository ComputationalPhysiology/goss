// Copyright (C) 2006-2012 Ola Skavhaug
//
// This file is part of GOSS.
//
// GOSS is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GOSS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GOSS. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Johan Hake 2012

#include <stdio.h>

#include "log.h"
#include "ESDIRK4O32.h"

using namespace goss;

//-----------------------------------------------------------------------------
ESDIRK4O32::ESDIRK4O32()
  : AdaptiveImplicitSolver(),
    gamma(0.43586652150845899942),
    a21(gamma),
    a22(gamma),
    a31((-4*gamma*gamma+6*gamma-1)/(4*gamma)),
    a32((-2*gamma+1)/(4*gamma)),
    a33(gamma),
    a41((6*gamma-1)/(12*gamma)),
    a42(-1/((24*gamma-12)*gamma)),
    a43((-6*gamma*gamma+6*gamma-1)/(6*gamma-3)),
    a44(gamma),
    b1(a41),
    b2(a42),
    b3(a43),
    b4(a44),
    bh1(a31),
    bh2(a32),
    bh3(a33),
    c2(2.0*gamma),
    c3(1.0),
    c4(1.0),
    z1(0), z2(0), z3(0), z4(0), yn(0), yh(0)
{
  _iord = 3;
}
//-----------------------------------------------------------------------------
ESDIRK4O32::ESDIRK4O32(std::shared_ptr<ODE> ode)
  : AdaptiveImplicitSolver(),
    gamma(0.43586652150845899942),
    a21(gamma),
    a22(gamma),
    a31((-4*gamma*gamma+6*gamma-1)/(4*gamma)),
    a32((-2*gamma+1)/(4*gamma)),
    a33(gamma),
    a41((6*gamma-1)/(12*gamma)),
    a42(-1/((24*gamma-12)*gamma)),
    a43((-6*gamma*gamma+6*gamma-1)/(6*gamma-3)),
    a44(gamma),
    b1(a41),
    b2(a42),
    b3(a43),
    b4(a44),
    bh1(a31),
    bh2(a32),
    bh3(a33),
    c2(2.0*gamma),
    c3(1.0),
    c4(1.0),
    z1(0), z2(0), z3(0), z4(0), yn(0), yh(0)
{
  _iord = 3;
  attach(ode);
}
//-----------------------------------------------------------------------------
ESDIRK4O32::ESDIRK4O32(const ESDIRK4O32& solver)
  : AdaptiveImplicitSolver(solver),
    gamma(0.43586652150845899942),
    a21(gamma),
    a22(gamma),
    a31((-4*gamma*gamma+6*gamma-1)/(4*gamma)),
    a32((-2*gamma+1)/(4*gamma)),
    a33(gamma),
    a41((6*gamma-1)/(12*gamma)),
    a42(-1/((24*gamma-12)*gamma)),
    a43((-6*gamma*gamma+6*gamma-1)/(6*gamma-3)),
    a44(gamma),
    b1(a41),
    b2(a42),
    b3(a43),
    b4(a44),
    bh1(a31),
    bh2(a32),
    bh3(a33),
    c2(2.0*gamma),
    c3(1.0),
    c4(1.0),
    z1(solver.num_states()), z2(solver.num_states()),
    z3(solver.num_states()), z4(solver.num_states()),
    yn(solver.num_states()), yh(solver.num_states())
{
  _iord = 3;
}
//-----------------------------------------------------------------------------
ESDIRK4O32::~ESDIRK4O32()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void  ESDIRK4O32::attach(std::shared_ptr<ODE> ode)
{

  // Use base classes to actually attach ode
  ImplicitODESolver::attach(ode);

  if (ode->is_dae())
    goss_error("ESDIRK4o32.cpp",
	       "attaching ode",
	       "cannot integrate a DAE ode with an explicit solver.");

  z1.resize(num_states());
  z2.resize(num_states());
  z3.resize(num_states());
  z4.resize(num_states());
  yn.resize(num_states());
  yh.resize(num_states());

}
//-----------------------------------------------------------------------------
void ESDIRK4O32::reset()
{
  // Reset counters
  nfevals = 0;
  ndtsa = 0;
  ndtsr = 0;
  _ldt = -1;

  // Reset bases
  AdaptiveImplicitSolver::reset();
}
//-----------------------------------------------------------------------------
void ESDIRK4O32::forward(double* y, double t, double interval)
{
  // NB sjekk definisjonen av prev vs hvilken dt som sendes til NewtonSolve!

  _t = t;
  uint i;
  const double t_end = _t + interval;

  double* ret_ptr = y;
  double* swap, *yn0;

  // Use the raw pointer enabling pointer swap
  yn0 = &yn[0];

  //if (ldt>0)
  //{
  //    // Calculates a local time step that is <= ldt,
  //    // and that divides dt into equally sized steps:
  //
  //    n_steps = (int) ceil(dt/ldt -1E-12);
  //    dt /= n_steps; //could have used new variable, but dt was available...
  //}

  _dt = _ldt;
  bool step_ok; // done = false;
  reached_tend = false;

  if (_dt < 0)
  {
    _ode->eval(y, _t, &z1[0]);
    _dt = dtinit(_t, y, yn0, &z1[0], &z2[0], _iord);
    nfevals += 2;
  }

#ifdef DEBUG
  // log data
  dt_v.push_back(_dt);
#endif

  while (!reached_tend)
  {

    // Computes the first node explicitly
    _ode->eval(y, _t, &z1[0]);
    nfevals += 1;

    // printf("z2=");
    // Computes the second node implicitly
    for (i = 0; i < num_states(); ++i)
    {
      _prev[i] = a21*z1[i];
      z2[i] = 0.0; //z1[i]-y[i]; // Ola comment: What happens here??
      //printf("%1.2e, ",z2[i]);
    }

    step_ok = newton_solve(&z2[0], &_prev[0], &y[0], _t + c2*_dt, _dt, a22, \
			   always_recompute_jacobian);

    // Need to check if the newton solver is converged.
    // If not, we half the stepsize and try again
    if (!step_ok)
    {
      _recompute_jacobian = true;
      _dt /= 2.0;
#ifdef DEBUG
      log_data(_dt, false);
      printf("Second node failed, new dt = %1.6e\n", _dt);
#endif
      continue;
    }

    //printf("z2=");
    for (i = 0; i < num_states(); ++i)
      z3[i] = z2[i] + y[i];

    //printf(" before eval\n");
    _ode->eval(&z2[0], _t + c2*_dt, &_f1[0]);
    nfevals += 1;

    // Computes the third node, implicitly
    // Use pointer swap instead of copy!!
    //printf("z2=");
    for (i = 0; i < num_states(); ++i)
    {
      //printf("%1.2e, ",f1[i]);
      z2[i] = _f1[i];
      _prev[i] = a31*z1[i] + a32*z2[i];
    }
    //printf(" after eval\n");

    step_ok = newton_solve(&z3[0], &_prev[0], y, _t + c3*_dt, _dt, a33,
			   always_recompute_jacobian);

    // Need to check if the newton solver is converged.
    // If not, we half the stepsize and try again
    if (!step_ok)
    {
      _recompute_jacobian=true;
      _dt /= 2.0;
#ifdef DEBUG
      log_data(_dt, false);
      printf("Third node failed, new dt = %1.6e\n", _dt);
#endif
      continue;
    }

    for (i = 0; i < num_states(); ++i)
      z3[i] += y[i];

    _ode->eval(&z3[0], _t + c3*_dt, &_f1[0]);
    nfevals += 1;

    // Computes the error estimate
    // Use pointer swap instead of copy!!
    for (i = 0; i < num_states(); ++i)
    {
      z3[i] = _f1[i];
      yh[i] = y[i] + _dt*(bh1*z1[i] + bh2*z2[i] + bh3*z3[i]);
    }

    // Computes the fourth node, implicitly
    for (i = 0; i < num_states(); ++i)
    {
      z4[i] = yh[i] - y[i];
      _prev[i] = a41*z1[i] + a42*z2[i] + a43*z3[i];
    }

    step_ok = newton_solve(&z4[0], &_prev[0], y, _t + c4*_dt, _dt, a44,
			   always_recompute_jacobian);

    // Need to check if the newton solver is converged.
    // If not, we half the stepsize and try again
    if (!step_ok)
    {
      _recompute_jacobian=true;
      _dt /= 2.0;
#ifdef DEBUG
      log_data(_dt, false);
      printf("Fourth node failed, new dt = %1.6e\n", _dt);
#endif
      continue;
    }

    for (i = 0; i < num_states(); ++i)
      z4[i] += y[i];

    _ode->eval(&z4[0], _t + c4*_dt, &_f1[0]);
    nfevals += 1;

    for (i = 0; i < num_states(); ++i)
    {
      //printf("i=%d, dt=%1.4e, b1=%1.4e, b2=%1.4e, b3=%1.4e, b4=%1.4e\n",i,dt,b1,b2,b3,b4);
      //printf("yn=%1.4e, y=%1.4e, z1=%1.4e, z2=%1.4e, z3=%1.4e, z4=%1.4e\n",yn[i],y[i],z1[i] ,z2[i] ,z3[i] ,f1[i]);
      yn[i] = y[i] + _dt*(b1*z1[i] + b2*z2[i] + b3*z3[i] + b4*_f1[i]);
      yh[i] -= yn[i];//use this as the error vector
    }

    new_time_step(y, yn0, &yh[0], t_end);//,blog);
    //recompute_jacobian=true;

    //done=blog[1];
#ifdef DEBUG
    log_data(_dt, step_accepted);
#endif
    if (step_accepted)
    {
      swap = y;
      y = yn0;
      yn0 = swap;
#ifdef DEBUG
      if (single_step_mode)
      {
        if (ret_ptr != y)
	{
          for (i = 0; i < num_states(); ++i)
            ret_ptr[i] = y[i];

          yn0 = y;
        }

        return;
      }
#endif
    }
  }

  // This is a copy to ensure that the input Ptr contains the final solution
  // This can probably be done in a more elegant way
  if (ret_ptr != y)
  {
    for (i = 0; i < num_states(); ++i)
      ret_ptr[i] = y[i];

    yn0 = y;
    //printf("Accepted at t=%1.4e, with dt=%1.4e\n",t,dt);
  }

#ifdef DEBUG
  dt_v.pop_back();
  //cout << "accepted,rejected="<<log[0]<<","<<log[1]<<endl;
  printf("ESDIRK4O32 done with comp_jac = %d and rejected = %d\n", _jac_comp, _rejects);
#endif

  //printf("Jumped out of the time loop at t = %1.6e\n",t);
}
//-----------------------------------------------------------------------------
