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

#include "gtest/gtest.h"

using namespace goss;


// Test fixture class template.
template <class O, class S>
class ODESolverTest : public testing::Test {
protected:
  ODESolverTest() : solver(new S(new O())) {}

  virtual ~ODESolverTest() { delete solver; }

  // ODE and solver
  ODESolver* solver;
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
class RLTester : public ODESolverTest<VDP, S>
{
public :
  RLTester() : ODESolverTest<VDP, S>() {}
};

// The list of ODEs we want to test.
typedef testing::Types<Arenstorf, Brusselator, NonLinOscillator, EulerRigidBody, \
		       FG, Robertson, Sin, VDP, SaltzLorenz> ODEs;

// Different list of Solvers
typedef testing::Types<ExplicitEuler, RK2, RK4, RKF32> ExplicitODESolvers; 
typedef testing::Types<ImplicitEuler, ESDIRK4O32> ImplicitODESolvers;
typedef testing::Types<RL, GRL1, GRL2> RLODESolvers;

TYPED_TEST_CASE(ODETester, ODEs);
TYPED_TEST_CASE(ExplicitTester, ExplicitODESolvers);
TYPED_TEST_CASE(ImplicitTester, ImplicitODESolvers);
TYPED_TEST_CASE(RLTester, RLODESolvers);

// Run all included 
TYPED_TEST(ODETester, IntegrationTest) 
{

  // Run coarse simulation
  this->run_ode(0.0001, 10.0, this->x_coarse);

  // Run fine simulation
  this->run_ode(0.00001, 10.0, this->x_fine);

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
