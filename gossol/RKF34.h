#ifndef RKF34_h_IS_INCLUDED
#define RKF34_h_IS_INCLUDED

#include "AdaptiveExplicitSolver.h"
#include <math.h>
#include <cstdlib>

namespace gossol 
{

  class RKF34 : public AdaptiveExplicitSolver {
    public:
      long nfevals, ndtsa, ndtsr; // A counter for the nunber of right hand side evaluations (nfevals) and the number of accepted and rejected timesteps (ndtsa, ndtsr)

      RKF34() {};
      RKF34 (gossol::ODE* ode_, double _ldt=-1.0) { 
        attach(ode_);
        ldt = _ldt; 
        init();
        AdaptiveExplicitSolver::init();
      } 

      RKF34(double _ldt) {
        ldt = _ldt;
        init();
        AdaptiveExplicitSolver::init();
      }

      void init(){
        iord    = 3;
      }

      virtual void attach(gossol::ODE* ode_);
      void forward(double* y, double t, double interval);

      void logData(double dt, bool accepted);
      void dtVector(gossol::DoubleVector *res);
      void acceptedVector(gossol::DoubleVector *res);

      ~RKF34 (){
        swap = NULL;
        free(ki);
        free(f1);
        free(f2f5);
        free(f3);
        free(f4);
        free(e);
        free(yn);
      };

    protected: 
      double a21,a31,a32,a41,a42,a43,a51,a53,a54;// RK coefficients
      double b1, b3, b4; // RK weights
      double bh1,bh3,bh4,bh5; // RK weights
      double d1, d3, d4, d5;// Error weights
      double c1,c2,c3,c4; // RK nodes
      long nbytes; // System size in bytes
      double *ki,*f1,*f2f5,*f3,*f4, *yn, *e, *swap, *retPtr;// state derivative, allocated in attach(ode)
      int itol;//parameter for scalar or vector tolerance computing
  };

}
#endif
