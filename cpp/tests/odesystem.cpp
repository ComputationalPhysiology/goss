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

#include <string>
#include <cmath>
#include <memory>
#include <boost/scoped_ptr.hpp>

#include "goss/goss.h"
#include "Winslow.h"
#include "WinslowCSE.h"
#include "WinslowCSEArray.h"
#include "Panfilov.h"
#include "PanfilovCSE.h"
#include "PanfilovCSEArray.h"

#include "gtest/gtest.h"

using namespace goss;

// Test fixture class template.
template <class O, class S>
class ODESystemSolverTest : public testing::Test {
protected:
  ODESystemSolverTest() : ode(new O()), solver(new S()){}

  virtual ~ODESystemSolverTest() {}

  // ODE and solver
  std::shared_ptr<ParameterizedODE> ode;
  std::shared_ptr<ODESolver> solver;
  boost::scoped_ptr<ODESystemSolver> system_solver;
  DoubleVector x_coarse, x_fine;

  double run_system(double dt, double tstop, DoubleVector& x, uint num_threads=0)
  {

    // Init ODESystemSolver
    system_solver.reset(new ODESystemSolver(x.n, solver, ode));

    // Reseting system solver
    system_solver->reset_default();

    // Set threads
    system_solver->set_num_threads(num_threads);

    // Fill solution vector with inital values
    system_solver->get_field_states(x.data.get());

    // Get init value
    const double init_value = x.data.get()[0];

    // Step solver and update field solution
    const uint nstep = std::ceil(tstop/dt - 1.0E-12);
    double t = 0.0;
    for (uint i = 0; i < nstep; i++)
    {
      system_solver->set_field_states(x.data.get());
      system_solver->forward(t, dt);
      system_solver->get_field_states(x.data.get());
      t += dt;
    }

    // Return init value for check
    return init_value;
  }
};

// Inherit to partially specialize
template<class O>
class ODETester : public ODESystemSolverTest<O, RK4>
{
public :
  ODETester() : ODESystemSolverTest<O, RK4>() {}
};

template<class S>
class ImplicitTester : public ODESystemSolverTest<Winslow, S>
{
public :
  ImplicitTester() : ODESystemSolverTest<Winslow, S>() {}
};

template<class S>
class OpenMPTester : public ODESystemSolverTest<Winslow, S>
{
public :
  OpenMPTester() : ODESystemSolverTest<Winslow, S>() {}
};

template<class S>
class ExplicitTester : public ODESystemSolverTest<WinslowCSE, S>
{
public :
  ExplicitTester() : ODESystemSolverTest<WinslowCSE, S>() {}
};

//template<class S>
//class RLTester : public ODESolverTest<VDP, S>
//{
//public :
//  RLTester() : ODESolverTest<VDP, S>() {}
//};

// The list of ODEs we want to test.
typedef testing::Types<Winslow, WinslowCSE, WinslowCSEArray, Panfilov, PanfilovCSE, \
		       PanfilovCSEArray> ODEs;

// Different list of Solvers
typedef testing::Types<RK2, RK4, RKF32> ExplicitODESolvers;
//typedef testing::Types<ImplicitEuler, ESDIRK4O32> ImplicitODESolvers;
typedef testing::Types<ImplicitEuler> ImplicitODESolvers;
//typedef testing::Types<RL, GRL1, GRL2> RLODESolvers;
//typedef testing::Types<Winslow> ParameterizedODEs;

TYPED_TEST_SUITE(ODETester, ODEs);
TYPED_TEST_SUITE(ExplicitTester, ExplicitODESolvers);
TYPED_TEST_SUITE(ImplicitTester, ImplicitODESolvers);
TYPED_TEST_SUITE(OpenMPTester, ImplicitODESolvers);

//TYPED_TEST_SUITE(RLTester, RLODESolvers);
//TYPED_TEST_SUITE(ParameterizedODETester, ParameterizedODEs);

// Run all included
TYPED_TEST(ODETester, IntegrationTest)
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  const double init_coarse = this->run_system(0.01, 10.0, this->x_coarse);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  const double init_fine = this->run_system(0.001, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

  // Check all of the data points are the same but not the same as the init
  for (uint i=1; i < this->x_fine.n; i++)
  {
    ASSERT_FLOAT_EQ(this->x_fine.data[0], this->x_fine.data[i]);
    ASSERT_FLOAT_EQ(this->x_coarse.data[0], this->x_coarse.data[i]);

    ASSERT_NE(init_fine, this->x_fine.data[i]);
    ASSERT_NE(init_coarse, this->x_coarse.data[i]);
  }
}

TYPED_TEST(ExplicitTester, ExplicitSolverTest)
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  const double init_coarse = this->run_system(0.01, 10.0, this->x_coarse);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  const double init_fine = this->run_system(0.001, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

  // Check all of the data points are the same but not the same as the init
  for (uint i=1; i < this->x_fine.n; i++)
  {
    ASSERT_FLOAT_EQ(this->x_fine.data[0], this->x_fine.data[i]);
    ASSERT_FLOAT_EQ(this->x_coarse.data[0], this->x_coarse.data[i]);

    ASSERT_NE(init_fine, this->x_fine.data[i]);
    ASSERT_NE(init_coarse, this->x_coarse.data[i]);
  }
}

TYPED_TEST(ImplicitTester, ImplicitSolverTest)
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  const double init_coarse = this->run_system(0.1, 10.0, this->x_coarse);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  const double init_fine = this->run_system(0.01, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

  // Check all of the data points are the same but not the same as the init
  for (uint i=1; i < this->x_fine.n; i++)
  {
    ASSERT_NEAR(this->x_fine.data[0], this->x_fine.data[i], 0.001);
    ASSERT_NEAR(this->x_coarse.data[0], this->x_coarse.data[i], 0.001);

    ASSERT_NE(init_fine, this->x_fine.data[i]);
    ASSERT_NE(init_coarse, this->x_coarse.data[i]);
  }
}

TYPED_TEST(OpenMPTester, OpenMPSolverTester)
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  const double init_coarse = this->run_system(0.1, 10.0, this->x_coarse, 4);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  const double init_fine = this->run_system(0.01, 10.0, this->x_fine, 4);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

  // Check all of the data points are the same but not the same as the init
  for (uint i=1; i < this->x_fine.n; i++)
  {
    ASSERT_NEAR(this->x_fine.data[0], this->x_fine.data[i], 0.001);
    ASSERT_NEAR(this->x_coarse.data[0], this->x_coarse.data[i], 0.001);

    ASSERT_NE(init_fine, this->x_fine.data[i]);
    ASSERT_NE(init_coarse, this->x_coarse.data[i]);
  }
}

//TYPED_TEST(RLTester, RLSolverTest)
//{
//
//  // Run coarse simulation
//  this->run_system(0.001, 10.0, this->x_coarse);
//
//  // Run fine simulation
//  this->run_system(0.0001, 10.0, this->x_fine);
//
//  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
//
//}
