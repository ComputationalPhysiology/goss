#ifndef ESDIRK4O32_h_IS_INCLUDED
#define ESDIRK4O32_h_IS_INCLUDED

#include "AdaptiveImplicitSolver.h"
#include <iostream>


namespace gossol 
{

  class ESDIRK4O32: public AdaptiveImplicitSolver{
    protected:

      double gamma;
      double a21, a22, a31, a32, a33, a41, a42, a43, a44;
      double b1, b2, b3, b4, bh1, bh2, bh3;
      double c2, c3, c4;  // RK coefficients
      double *z1, *z2, *z3, *z4, *yn, *yh, *swap, *retPtr; 


    public:
      long nfevals, ndtsa, ndtsr; // A counter for the nunber of right hand side evaluations (nfevals) and the number of accepted and rejected timesteps (ndtsa, ndtsr)

      ESDIRK4O32() {};
      ESDIRK4O32 (gossol::ODE* ode_, double _ldt=-1.0) { 
        attach(ode_);
        ldt = _ldt; 
        init();
        AdaptiveImplicitSolver::init();
      } 

      ESDIRK4O32(double _ldt) {
        ldt = _ldt;
        init();
        AdaptiveImplicitSolver::init();
      }

      void init(){
        iord    = 3;
      }


      virtual void attach(gossol::ODE* ode_);
      void forward(double* y, double t, double dt);

      ~ESDIRK4O32 () {
        swap = NULL;
        retPtr=NULL;
        free(z1);
        free(z2);
        free(z3);
        free(z4);
        free(yn);
        free(yh);

      }
  };

}
#endif
