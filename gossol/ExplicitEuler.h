#ifndef ExplicitEuler_h_IS_INCLUDED
#define ExplicitEuler_h_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>
#include <cstdlib>

namespace gossol 
{

  class ExplicitEuler : public ODESolver
  {
    public:
      ExplicitEuler(): ODESolver(0.0, 0.0), dFdt(0) {};
      ExplicitEuler (gossol::ODE* ode_, double _ldt=-1.0);
      ExplicitEuler(double _ldt);

      virtual void attach(gossol::ODE* ode_);
      void forward(double* y, double t, double interval);
      ~ExplicitEuler () { free(dFdt); }

    protected: 
      double* dFdt; // state derivative, allocated in attach(ode)
  };

}
#endif
