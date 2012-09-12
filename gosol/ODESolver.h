#ifndef ODESolver_h_IS_INCLUDED
#define ODESolver_h_IS_INCLUDED

#include "../ODE.h"
#include <iostream>

namespace gosol 
{

  class ODESolver {
    protected:
      double ldt; // local time step.
      double  dt; // variable local time step.
      gosol::ODE* ode;
      long num_tsteps;

    public:
      virtual ~ODESolver () { /* Do nothing */ }
      ODESolver (gosol::ODE* ode_, double ldt=0.0, double dt=0.0);
      ODESolver (double ldt=0.0, double dt=0.0);

      virtual void attach(gosol::ODE* ode_);
      virtual void forward(double* y, double t, double DT) = 0;
      void probe(double *y) const {ode->probe(y); }
      inline int ode_size() const { return ode->size(); }
      inline const gosol::ODE*  odesystem() const { return ode; }
      inline const double internalTimeStep() const { return ldt; }
      inline void setInternalTimeStep(double idt) { ldt=idt; }
  };

}
#endif
