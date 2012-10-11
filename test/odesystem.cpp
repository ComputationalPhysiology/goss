#include <string>
#include <cmath>

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
  ODESystemSolverTest() {}

  virtual ~ODESystemSolverTest() {}

  // ODE and solver
  ParameterizedODE* ode;
  ODESolver* solver;
  ODESystemSolver* system_solver;
  DoubleVector x_coarse, x_fine;
  
  void run_system(double dt, double tstop, DoubleVector& x, uint num_threads=0)
  {

    // Init ODESystemSolver
    ode = new O();
    solver = new S();
    system_solver = new ODESystemSolver(x.n, solver, ode);

    // Reseting system solver
    system_solver->reset_default();

    // Set threads
    system_solver->set_num_threads(num_threads);

    // Fill solution vector with inital values
    system_solver->get_field_states(x.data.get());

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

    // Clean up
    delete system_solver;
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

TYPED_TEST_CASE(ODETester, ODEs);
TYPED_TEST_CASE(ExplicitTester, ExplicitODESolvers);
TYPED_TEST_CASE(ImplicitTester, ImplicitODESolvers);
TYPED_TEST_CASE(OpenMPTester, ImplicitODESolvers);

//TYPED_TEST_CASE(RLTester, RLODESolvers);
//TYPED_TEST_CASE(ParameterizedODETester, ParameterizedODEs);

// Run all included 
TYPED_TEST(ODETester, IntegrationTest) 
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  this->run_system(0.01, 10.0, this->x_coarse);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  this->run_system(0.001, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
  
}

TYPED_TEST(ExplicitTester, ExplicitSolverTest) 
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  this->run_system(0.01, 10.0, this->x_coarse);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  this->run_system(0.001, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
}

TYPED_TEST(ImplicitTester, ImplicitSolverTest) 
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  this->run_system(0.1, 10.0, this->x_coarse);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  this->run_system(0.01, 10.0, this->x_fine);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
  
}

TYPED_TEST(OpenMPTester, OpenMPSolverTester) 
{

  // Run coarse simulation
  const uint nodes(100);
  this->x_coarse.data.reset(new double[nodes]);
  this->x_coarse.n = nodes;
  this->run_system(0.1, 10.0, this->x_coarse, 4);

  // Run fine simulation
  this->x_fine.data.reset(new double[nodes]);
  this->x_fine.n = nodes;
  this->run_system(0.01, 10.0, this->x_fine, 4);

  ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
  
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
