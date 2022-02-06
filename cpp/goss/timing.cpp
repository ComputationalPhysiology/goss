// Copyright (C) 2003-2011 Anders Logg
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
// First added:  2003-12-21
// Last changed: 2014-04-07

// Uncomment this for testing std::clock
//#define _WIN32

#ifdef _WIN32
#include <boost/timer.hpp>
#else
#include <sys/time.h>
#endif

#include "log.h"
#include "LogManager.h"
#include "timing.h"

// We use std::timer (std::clock) on Windows and otherwise the
// platform-dependent (but higher-precision) gettimeofday from
// <sys/time.h>. Note that in the latter case, the timer is not
// reset to zero at the start of the program so time() will not
// report total CPU time, only the difference makes sense.

namespace goss
{
#ifdef _WIN32
  std::timer __global_timer;
  std::timer __tic_timer;
#else
  double __tic_timer;
#endif
}

using namespace goss;

//-----------------------------------------------------------------------
void goss::tic()
{
#ifdef _WIN32
  goss::__tic_timer.restart();
#else
  goss::__tic_timer = time();
#endif
}
//-----------------------------------------------------------------------------
double goss::toc()
{
#ifdef _WIN32
  return __tic_timer.elapsed();
#else
  return time() - __tic_timer;
#endif
}
//-----------------------------------------------------------------------------
double goss::time()
{
#ifdef _WIN32
  return goss::__global_timer.elapsed();
#else
  struct timeval tv;
  struct timezone tz;
  if (gettimeofday(&tv, &tz) != 0)
  {
    goss_error("timing.cpp",
	       "return current time",
	       "Call to gettimeofday() failed");
  }
  return static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec)*1e-6;
#endif
}
//-----------------------------------------------------------------------------
void goss::list_timings(bool reset)
{
  // Optimization
  if (!LogManager::logger.is_active())
    return;

  LogManager::logger.list_timings(reset);
}
//-----------------------------------------------------------------------------
Table goss::timings(bool reset)
{
  return LogManager::logger.timings(reset);
}
//-----------------------------------------------------------------------------
double goss::timing(std::string task, bool reset)
{
  return LogManager::logger.timing(task, reset);
}
//-----------------------------------------------------------------------------
