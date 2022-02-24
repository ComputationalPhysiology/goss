/* -*- C -*- */
// Copyright (C) 2012 Johan Hake
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

#ifndef __GOSS_H
#define __GOSS_H

// GOSS interface

// Common
#include <goss/types.h>
#include <goss/constants.h>
#include <goss/DoubleVector.h>

// log
#include <goss/log.h>
#include <goss/LogStream.h>
#include <goss/Progress.h>
#include <goss/Table.h>
#include <goss/LogLevel.h>

// Timing
#include <goss/timing.h>
#include <goss/Timer.h>

// ODEs
#include <goss/ODE.h>
#include <goss/ParameterizedODE.h>

// ODESolvers
#include <goss/ODESolver.h>
#include <goss/ImplicitODESolver.h>
#include <goss/AdaptiveExplicitSolver.h>
#include <goss/AdaptiveImplicitSolver.h>
#include <goss/ExplicitEuler.h>
#include <goss/RK4.h>
#include <goss/RK2.h>
//#include <goss/RKF34.h>
#include <goss/RKF32.h>
#include <goss/BasicImplicitEuler.h>
#include <goss/ImplicitEuler.h>
#include <goss/ThetaSolver.h>
//#include <goss/ImplicitMidPointSolver.h>
//#include <goss/SDIRK2O2.h>
#include <goss/ESDIRK4O32.h>
#include <goss/ESDIRK23a.h>
#include <goss/GRL1.h>
#include <goss/GRL2.h>
#include <goss/RL1.h>
#include <goss/RL2.h>

// ODESystem
#include <goss/ODESystemSolver.h>

#endif
