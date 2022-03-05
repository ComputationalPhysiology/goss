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

#include "gtest/gtest.h"

#include <string>
#include <cmath>

#include "goss/goss.h"
#include "NonLinOscillator.h"
#include "Arenstorf.h"
#include "Brusselator.h"
#include "EulerRigidBody.h"
#include "FG.h"
#include "Robertson.h"
#include "SaltzLorenz.h"
#include "Sin.h"
#include "VDP.h"
#include "van_der_Pol.h"
#include "Winslow.h"
#include "WinslowNoIntermediates.h"
#include "WinslowCSE.h"
#include "WinslowCSEArray.h"
#include "Panfilov.h"
#include "PanfilovNoIntermediates.h"
#include "PanfilovCSE.h"
#include "PanfilovCSEArray.h"
#include "fhn.h"
#include "beeler_reuter_1977.h"

using namespace goss;


// Test fixture class template.
template <class O, class S>
class ODESolverTest : public testing::Test {
protected:
  ODESolverTest() : ode(new O()), solver(new S(ode)) {}

  virtual ~ODESolverTest() { }

  // ODE and solver
  std::shared_ptr<ODE> ode;
  std::shared_ptr<ODESolver> solver;
  DoubleVector x_coarse, x_fine;

  void run_ode(double dt, double tstop, DoubleVector& x)
  {

    // Reseting solver
    solver->reset();
    solver->get_ode()->get_ic(&x);

    const uint nstep = std::ceil(tstop/dt - 1.0E-12);

    double t = 0.0;

    for (uint i = 0; i < nstep; i++)
    {
      solver->forward(x.data.get(), t, dt);
      t += dt;
      /*
      if (std::fmod(t, .1) < dt)
	info("t %.2f: y: %.2f, z: %.2f", t, x.data.get()[0], x.data.get()[1]);
      */
    }
  }
};

// Inherit to partially specialize
template<class O>
class ODETester : public ODESolverTest<O, RK4>
{
public :
  ODETester() : ODESolverTest<O, RK4>() {}
};

template<class O>
class DAETester : public ODESolverTest<O, ThetaSolver>
{
public :
  DAETester() : ODESolverTest<O, ThetaSolver>() {}
};

template<class S>
class ImplicitTester : public ODESolverTest<EulerRigidBody, S>
{
public :
  ImplicitTester() : ODESolverTest<EulerRigidBody, S>() {}
};

template<class S>
class ExplicitTester : public ODESolverTest<Arenstorf, S>
{
public :
  ExplicitTester() : ODESolverTest<Arenstorf, S>() {}
};

template<class S>
class RLTester : public ODESolverTest<Fhn, S>
{
public :
  RLTester() : ODESolverTest<Fhn, S>() {}
};

template<class S>
class RLTester2 : public ODESolverTest<Beeler_reuter_1977, S>
{
public :
  RLTester2() : ODESolverTest<Beeler_reuter_1977, S>() {}
};

template<class O>
class ParameterizedODETester : public ODESolverTest<O, RKF32>
{
public :
  ParameterizedODETester() : ODESolverTest<O, RKF32>() {}
};

// The list of ODEs we want to test.
typedef testing::Types<Arenstorf, Brusselator, NonLinOscillator, EulerRigidBody, \
		       FG, Robertson, Sin, VDP, SaltzLorenz> ODEs;

typedef testing::Types<VanDerPol> DAEs;

// Different list of Solvers
typedef testing::Types<ExplicitEuler, RK2, RK4> ExplicitODESolvers;
//typedef testing::Types<ImplicitEuler, SDIRK> ImplicitODESolvers;
typedef testing::Types<BasicImplicitEuler, ESDIRK23a, ImplicitEuler, ThetaSolver> ImplicitODESolvers;
typedef testing::Types<RL1, GRL1, RL2, GRL2> RLODESolvers;
typedef testing::Types<Winslow, WinslowNoIntermediates, WinslowCSE,
		       WinslowCSEArray, Panfilov, PanfilovNoIntermediates,
		       PanfilovCSE, PanfilovCSEArray> ParameterizedODEs;

TYPED_TEST_SUITE(ODETester, ODEs);
TYPED_TEST_SUITE(DAETester, DAEs);
TYPED_TEST_SUITE(ExplicitTester, ExplicitODESolvers);
TYPED_TEST_SUITE(ImplicitTester, ImplicitODESolvers);
TYPED_TEST_SUITE(RLTester, RLODESolvers);
TYPED_TEST_SUITE(ParameterizedODETester, ParameterizedODEs);

// Run all included
TYPED_TEST(ODETester, IntegrationTest)
{

  // Run coarse simulation
  this->run_ode(0.0001, 10.0, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.00001, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

}

TYPED_TEST(DAETester, DAESolverTest)
{

  // Run coarse simulation
  this->run_ode(0.01, 0.800, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.001, 0.800, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

}

TYPED_TEST(ExplicitTester, ExplicitSolverTest)
{

  // Run coarse simulation
  this->run_ode(0.0001, 10.0, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.00001, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
}

TYPED_TEST(ImplicitTester, ImplicitSolverTest)
{

  // Run coarse simulation
  this->run_ode(0.1, 10.0, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.01, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

}

TYPED_TEST(RLTester, RLSolverTest)
{

  // Run coarse simulation
  this->run_ode(0.001, 10.0, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.0001, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

}

TYPED_TEST(ParameterizedODETester, ParameterizedODETest)
{

  // Local version
  ParameterizedODE& lode = dynamic_cast<ParameterizedODE&>(*this->ode);

  // Run coarse simulation
  this->run_ode(0.1, 10.0, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.01, 10.0, this->x_fine);

  double old_fine = this->x_fine.data[0];
  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);

  // Change parameters

  std::string param = lode.get_parameter_names()[0];
  const double value = lode.get_parameter(param);
  lode.set_parameter(param, value*0.9);

  ASSERT_EQ(lode.get_parameter(param), value*0.9);

  // Set field parameters
  std::vector<double> values(lode.num_field_parameters(), 0.0);
  for (uint i = 0; i < lode.num_field_parameters(); i++)
  {
    param = lode.get_field_parameter_names()[i];
    values[i] = lode.get_parameter(param)*0.9;
  }

  // Set field parameters if any
  if (lode.num_field_parameters() > 0)
  {
    lode.set_field_parameters(&values[0]);
    ASSERT_EQ(lode.get_parameter(param), values[values.size()-1]);
  }

  // Run coarse simulation
  this->run_ode(0.1, 100.0, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.01, 100.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
  ASSERT_FALSE(std::fabs(old_fine - this->x_fine.data[0]) < .001);

}
