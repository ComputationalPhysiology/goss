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
template <class T>
class ODETest : public testing::Test {
protected:
  ODETest() : ode(new T()), solver(new ExplicitEuler()) {}

  virtual ~ODETest() { delete ode; delete solver;}

  // ODE and solver
  ODE* ode;
  ODESolver* solver;
  DoubleVector x_coarse, x_fine;
  
  void run_ode(const double dt, DoubleVector& x)
  {
    ode->getIC(&x);
    const double tstop = .1;
    const int nstep = ceil(tstop/dt - 1.0E-12); 
    double t = 0.0;
    for (int i = 0; i < ode->size(); i++)
      printf("var %d: %f", i, x.data[i]);
    
    for (int i = 0; i < nstep; i++)
    {
      solver->forward(x.data, t, dt);
      t += dt;
    }
  }
};

// The list of types we want to test.
typedef testing::Types<NonLinOscillator, Arenstorf, Brusselator, EulerRigidBody, FG, Robertson, SaltzLorenz, Sin, VDP> ODEs;

TYPED_TEST_CASE(ODETest, ODEs);

TYPED_TEST(ODETest, IntegrationTest) {

  // Run coarse simulation
  this->run_ode(0.05, this->x_coarse);

  // Run fine simulation
  //this->run_ode(0.001, this->x_fine);

  //ASSERT_NEAR(this->x_fine.data[0], this->x_coarse.data[0], 1.0);
  
}