#ifndef ODESolver_h_IS_INCLUDED
#define ODESolver_h_IS_INCLUDED

#include <iostream>

#include "ODE.h"

namespace goss 
{

  class ODESolver {
    protected:
      double ldt; // local time step.
      double  dt; // variable local time step.
      goss::ODE* ode;
      long num_tsteps;

    public:
      virtual ~ODESolver () { /* Do nothing */ }
      ODESolver (goss::ODE* ode_, double ldt=0.0, double dt=0.0);
      ODESolver (double ldt=0.0, double dt=0.0);

      virtual void attach(goss::ODE* ode_);
      virtual void forward(double* y, double t, double DT) = 0;
      void probe(double *y) const {ode->probe(y); }
      inline int ode_size() const { return ode->size(); }
      inline const goss::ODE*  odesystem() const { return ode; }
      inline const double internalTimeStep() const { return ldt; }
      inline void setInternalTimeStep(double idt) { ldt=idt; }
  };

}
#endif
