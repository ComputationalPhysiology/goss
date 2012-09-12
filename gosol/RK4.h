#ifndef RK4_h_IS_INCLUDED
#define RK4_h_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>

namespace gosol 
{

  class RK4 : public ODESolver {
    public:
      RK4() {};
      RK4(double _ldt);
      RK4(gosol::ODE* ode_, double _ldt=-1.0);
      ~RK4();

      virtual void attach(gosol::ODE* ode_);
      void forward(double* y, double t, double interval);
      inline void add(double* x, double* y, double a, double* z);

    protected: 
      double *k1, *k2, *k3, *k4, *tmp;// state derivative, allocated in attach(ode)

  };

  inline void RK4:: add(double* x, double* y, double a, double* z)
  {
    for (int i=0; i < ode->size(); ++i) {
      x[i] = y[i];
      x[i] += a*z[i];
    }
  }

}
#endif
