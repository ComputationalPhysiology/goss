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

#include "AdaptiveExplicitSolver.h"
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace goss;

//-----------------------------------------------------------------------------
AdaptiveExplicitSolver::AdaptiveExplicitSolver()
  : ODESolver(), num_accepted(0), num_rejected(0), _t(0.), _ldt(0.1),
    _dt(0.1), _dt_prev(0.), _atol(1.e-5), _rtol(1.e-8), _iord(1),
    facmin(0.5), facmax(2.0), facmaxb(facmax), stabfac(0.9),
    step_accepted(false), reached_tend(false),
    _itol(0), dt_v(0), accept_v(0), single_step_mode(false)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
AdaptiveExplicitSolver::AdaptiveExplicitSolver(const AdaptiveExplicitSolver& solver)
  : ODESolver(solver), num_accepted(solver.num_accepted),
    num_rejected(solver.num_rejected), _t(solver._t), _ldt(solver._ldt), _dt(solver._dt),
    _dt_prev(solver._dt_prev),_atol(solver._atol), _rtol(solver._rtol),
    _iord(solver._iord), facmin(solver.facmin), facmax(solver.facmax),
    facmaxb(solver.facmaxb), stabfac(solver.stabfac),
    step_accepted(solver.step_accepted), reached_tend(solver.reached_tend),
    _itol(solver._itol), dt_v(solver.dt_v), accept_v(solver.accept_v),
    single_step_mode(solver.single_step_mode)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver::reset()
{
  single_step_mode = false;
  reached_tend     = false;
  num_accepted = 0;
  num_rejected = 0;
  _atol    = 1.0e-5;
  _rtol    = 1.0e-8;
  _itol    = 0;

  // We can not choose the next timestep more then
  // half of the previous timestep
  facmin  = 0.5;

  // We can not choose the next timestep more then
  // double of the previous timestep
  facmaxb = 2.0;
  facmax  = facmaxb;
  stabfac = 0.9; // std::pow(0.25,1/(_iord+1));
}
//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::dtinit(double t, double* y0, double* y1,
				      double* f0, double* f1, double iord)
{
  // Compute a first guess for explicit Euler as
  // H = 0.01 * norm(y0)/(norm(f0)
  // The increment for explicit Euler is small compared to the solution
  // We assume that y0 and f0 are computed.
  // y1 and f1 are just pointers to contigous memory which this
  // function borrows

  uint i;
  double dnf = 0.0;
  double dny = 0.0;
  double sk, dt, tmp;
  if (_itol == 0)
  {
    for (i = 0; i < num_states(); ++i)
    {
      sk   = _atol + _rtol*fabs(y0[i]);
      tmp  = f0[i]/sk;
      dnf += tmp*tmp;
      tmp  = y0[i]/sk;
      dny += tmp*tmp;
    }
  }
  else
  {
    std::cout << " Not implemented yet " << std::endl;
    //for (i=0; i<num_states();++i){
    //    sk   = _atol[i] + _rtol[i]*fabs(y0[i]);
    //    dnf += pow(f0[i]/sk,2);
    //    dny += pow(y0[i]/sk,2);
    //}
    //for i in xrange(0,n):
    //    sk   = _atol[i] + _rtol[i]*math.fabs(y0[i])
    //    dnf += (f0[i]/sk)**2
    //    dny += (y0[i]/sk)**2
    //print "dnf = ", dnf, " dny = ", dny
  }

  if (dnf <= 1.0e-10 || dny <= 1.0e-10)
    dt = 1.0e-6;
  else
    dt = 0.01*std::sqrt(dny/dnf);

  // Should we have a dt_max??
  // Perform an explicit Euler step
  for (i = 0; i < num_states(); ++i)
    y1[i] = y0[i] + dt*f0[i];

  _ode->eval(y1, t + dt, f1);

  // Estimate the second derivative of the solution
  double der2 = 0.0;
  if (_itol == 0)
  {
    for (i = 0; i < num_states(); ++i)
    {
      sk    = _atol + _rtol*fabs(y1[i]);
      tmp   = (f1[i] - f0[i])/sk;
      der2 += tmp*tmp;
    }
  }
  else
  {
    std::cout << " Not implemented yet" << std::endl;
    //for (i=0;i<num_states();++i){
    //    sk    = _atol[i] + _rtol[i]*fabs(y1[i]);
    //    der2 += pow(((f1[i]-f0[i])/sk),2);
    //}
    //for i in xrange(0,n):
    //     sk    = _atol[i] + _rtol[i]*math.fabs(y0[i])
    //     der2 += ((f1[i]-f0[i])/sk)**2
  }
  der2 = std::sqrt(der2)/dt;

  // Step size is computed such that
  // dt**iord*max(norm(f0), norm(der2)) = 0.01

  double der12;
  if (fabs(der2) >= std::sqrt(dnf))
    der12 = fabs(der2);
  else
    der12 = std::sqrt(dnf);

  double dt1;
  if (der12 <= 1.0e-15)
  {
    if (fabs(dt)*1.0e-3 < 1.0e-6)
      dt1 = 1.0e-6;
    else
      dt1 = fabs(dt)*1.0e-3;
  }
  else
    dt1 = std::pow((0.01/der12),(1.0/iord));

  if (100*fabs(dt)<dt1)
    return 100*fabs(dt);
  else
    return dt1;
}
//-----------------------------------------------------------------------------
void AdaptiveExplicitSolver::new_time_step(double* y, double* yn, double* e, double t_end)
{
  // FIXME: Force this method to return the new time step
  double err = 0.0;
  double sk;

  // A way to check if we are at t_end.
  double eps = 1e-14;
  double max;
  double yi_abs, yni_abs, tmp;
  _dt_prev = _dt;

  for (uint i = 0; i < num_states(); ++i)
  {
    yi_abs  = fabs(y[i]);
    yni_abs = fabs(yn[i]);
    max     = yi_abs > yni_abs ? yi_abs : yni_abs;
    sk      = _atol + _rtol*max;
    tmp     = e[i]/sk;
    err += tmp*tmp;
  }

  err = std::sqrt(err/num_states());

  // Sanity check to ensure algorithm does not break down.
  // If we encounter strange values in the error, set it large to enforce a
  // small time step for the next iteration.
  if (std::isnan(err) || std::isinf(err))
    err = 2000.0;

  // If the error is smaller than 1 the timestep is accepted, and we advance
  // If not, the timestep is rejected
  if (err <= 1.0)
  {
    _t += _dt;
    num_accepted += 1;
    step_accepted = true;

    if (fabs(_t - t_end) < eps)
      reached_tend = true;
  }
  else
  {
    num_rejected += 1;
    step_accepted = false;
  }

  // Computation of dtnew
  double fac = stabfac*std::pow(1.0/err, 1.0/(_iord + 1));
  if (facmin > fac)
    fac = facmin;

  if (fac > facmax)
    fac = facmax;

  // If the timestep is rejected, we prevent the next timestep from increasing
  if (!step_accepted)
    facmax = 1.0;
  else
    facmax = facmaxb;

  _dt *= fac;

  // If it is very likely that we will reach the end with two timesteps, we
  // set the first timestp to half the distance to avoid the last timestep
  // been very small
  if (_t + _dt >= t_end)
    _dt = t_end - _t;

  else if (_t + 1.5*_dt >= t_end)
    _dt = (t_end - _t)/2.0;

  else
    _ldt = _dt;// Saves the timestep to be used as initial guess for next macro step

  //std::cout << "dt = " << _dt << std::endl;
}
//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::get_current_time()
{
#ifndef DEBUG
  printf("NOT IN DEBUG MODE\n");
#endif
  return _t;
}
//-----------------------------------------------------------------------------
double AdaptiveExplicitSolver::get_current_time_step()
{
#ifndef DEBUG
  printf("NOT IN DEBUG MODE\n");
  return _dt;
#else
  return _dt_prev;
#endif
}
//-----------------------------------------------------------------------------
