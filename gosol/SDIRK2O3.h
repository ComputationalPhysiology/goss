#ifndef SDIRK2O3_h_IS_INCLUDED
#define SDIRK2O3_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <iostream>
#include <vector>

namespace gosol 
{

  class SDIRK2O3: public ImplicitODESolver{
    public:
      std::vector<int> newtonIter1;
      std::vector<int> newtonIter2;
      std::vector<int> newtonAccepted1;
      std::vector<int> newtonAccepted2;
      std::vector<double> dt_v;

      SDIRK2O3() {};
      SDIRK2O3 (gosol::ODE* ode_, double _ldt=-1.0) { 
        attach(ode_);
        ldt = _ldt; 
      } 

      SDIRK2O3(double _ldt) {
        ldt = _ldt;
      }

      virtual void attach(gosol::ODE* ode_);
      void forward(double* y, double t, double interval);
      ~SDIRK2O3 () {
        free(z1);
        free(z2);
      }

    protected:
      double ldt; // local time step.
      double gamma;
      double a11, a21, a22, b1, b2, c1, c2, d1, d2;// RK coefficients
      double *z1, *z2; //one intermediate stage g1
      bool justrefined;
  };

}
#endif
