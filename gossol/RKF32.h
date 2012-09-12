#ifndef RKF32_h_IS_INCLUDED
#define RKF32_h_IS_INCLUDED

#include "AdaptiveExplicitSolver.h"
#include <math.h>

namespace gossol 
{

  class RKF32 : public  AdaptiveExplicitSolver
  {
    public:
      long nfevals, ndtsa, ndtsr; // A counter for the nunber of right hand side evaluations (nfevals) and the number of accepted and rejected timesteps (ndtsa, ndtsr)

      RKF32() {};
      RKF32 (gossol::ODE* ode_, double _ldt=-1.0);
      RKF32(double _ldt);
      ~RKF32();

      void init();
      virtual void attach(gossol::ODE* ode_);
      void forward(double* y, double t_, double interval);

      void logData(double dt, bool accepted);
      void dtVector(gossol::DoubleVector *res);
      void acceptedVector(gossol::DoubleVector *res);

    protected: 
      double a21,a32;// RK coefficients
      double b1,b2,b3; // RK weights
      double bh1,bh2,bh3,bh4; // RK weights
      double d1,d2,d3,d4;//Error weights
      double c1,c2,c3; // RK nodes
      long nbytes; // System size in bytes
      double *ki,*k1,*k2,*k3,*k4, *yn, *e, *swap, *retPtr;// state derivative, allocated in attach(ode
      int itol;//parameter for scalar or vector tolerance computing
      bool first;
  };

}
#endif
