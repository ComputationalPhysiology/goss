#ifndef ImplicitEuler_h_IS_INCLUDED
#define ImplicitEuler_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <vector>

namespace goss 
{

  class ImplicitEuler: public ImplicitODESolver{
    public:
      std::vector<int> newtonIter1;
      std::vector<int> newtonAccepted1;
      std::vector<double> dt_v;

      ImplicitEuler() : z1(0) {};
      ImplicitEuler (goss::ODE* ode_, double _ldt=-1.0) : z1(0) { 
        attach(ode_);
        ldt = _ldt; 
      } 

      ImplicitEuler(double _ldt) : z1(0) {
        ldt = _ldt;
      }

      ~ImplicitEuler () {free(z1);}

      virtual void attach(goss::ODE* ode_);
      void forward(double* y, double t, double dt);

    protected:
      double ldt; // local time step.
      double* z1;
      bool justrefined;


  };

}
#endif
