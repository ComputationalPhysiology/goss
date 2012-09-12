#ifndef SDIRK2O2_h_IS_INCLUDED
#define SDIRK2O2_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <iostream>
#include <vector>


namespace gosol 
{

  class SDIRK2O2: public ImplicitODESolver
  {
    public:
      std::vector<int> newtonIter1;
      std::vector<int> newtonIter2;
      std::vector<int> newtonAccepted1;
      std::vector<int> newtonAccepted2;
      std::vector<double> dt_v;

      SDIRK2O2() {};
      SDIRK2O2 (gosol::ODE* ode_, double _ldt=-1.0);
      SDIRK2O2(double _ldt);
      ~SDIRK2O2();

      virtual void attach(gosol::ODE* ode_);
      void forward(double* y, double t, double interval);

    protected:
      double gamma;
      double a11, a21, a22, b1, b2, c1, c2;// RK coefficients
      double *z1, *z2; //one intermediate stage g1
      bool justrefined;
  };

}
#endif
