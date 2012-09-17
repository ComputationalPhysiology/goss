#ifndef ThetaSolver_h_IS_INCLUDED
#define ThetaSolver_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <vector>


namespace goss 
{

  class ThetaSolver: public ImplicitODESolver{
    public:
      double theta;
      std::vector<int> newtonIter1;
      std::vector<int> newtonAccepted1;
      std::vector<double> dt_v;

      ThetaSolver() {};
      ThetaSolver (goss::ODE* ode_, double _ldt=-1.0);
      ThetaSolver(double _ldt);
      ~ThetaSolver ();

      virtual void attach(goss::ODE* ode_);
      void forward(double* y, double t, double dt);

    protected:
      double ldt; // local time step.
      double *z1, *ft1;
      bool justrefined;

  };

}
#endif
