#ifndef RK2_H_IS_INCLUDED
#define RK2_H_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>

namespace gosol 
{

  class RK2 : public ODESolver {
    public:
      RK2() {};
      RK2(gosol::ODE* ode_);
      ~RK2();

      virtual void attach(gosol::ODE* ode_);
      void forward(double* y, double t, double dt);
      inline void add(double* x, double* y, double a, double* z);

    protected: 
      double *k1, *tmp;// state derivative, allocated in attach(ode)

  };

  inline void RK2:: add(double* x, double* y, double a, double* z)
  {
    for (int i=0; i < ode->size(); ++i) {
      x[i] = y[i];
      x[i] += a*z[i];
    }
  }

}
#endif
