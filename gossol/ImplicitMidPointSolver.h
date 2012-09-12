#ifndef ImplicitMidPointSolver_h_IS_INCLUDED
#define ImplicitMidPointSolver_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <vector>

namespace gossol 
{

  class ImplicitMidPointSolver: public ImplicitODESolver{
    public:
      double theta;
      std::vector<int> newtonIter1;
      std::vector<int> newtonAccepted1;
      std::vector<double> dt_v;

      ImplicitMidPointSolver() {};
      ImplicitMidPointSolver (ODE* ode_, double _ldt=-1.0) { 
        attach(ode_);
        ldt = _ldt; 
      } 

      ImplicitMidPointSolver(double _ldt) {
        ldt = _ldt;
      }

      ~ImplicitMidPointSolver () {free(z1);}

      virtual void attach(ODE* ode_);
      void forward(double* y, double t, double dt);

    protected:
      double ldt; // local time step.
      double *z1;//, *f1;
      bool justrefined;
  };

}
#endif
