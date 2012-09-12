#ifndef AdaptiveExplicitSolver_h_IS_INCLUDED
#define AdaptiveExplicitSolver_h_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>
#include <stdio.h>
#include <vector>

namespace gossol 
{

  class AdaptiveExplicitSolver : public ODESolver
  {
    public:
      AdaptiveExplicitSolver (gossol::ODE* ode_, double ldt=0.0, double dt=0.0);
      AdaptiveExplicitSolver(double ldt=0.0, double dt=0.0);
      virtual ~AdaptiveExplicitSolver (){};

      void init();
      /*virtual void attach(gossol::ODE* ode_);*/
      virtual void forward(double* y, double t, double interval) = 0;
      double getCurrentTime();
      double getCurrentTimeStep();
      long getNumAccepted(){return num_accepted;}
      long getNumRejected(){return num_rejected;}
      void setSingleStepMode(bool mode){single_step_mode=mode;}
      void setTol(double atol_, double rtol_=1.0e-8);
      void setIord(int iord_);
      double dtinit(double t, double* y0, double* y1, double* f0, double* f1, double iord);
      void newTimeStep(double* y, double* yn, double* e, double t_end);

    protected: 
      long num_accepted, num_rejected;//a log of 1) the numer of steps, 2)the number of rejected steps
      double t, dt_prev;
      double atol,  rtol, iord, facmin, facmax, facmaxb, stabfac; // local time step and tolerence.
      bool step_accepted, reached_tend;// a bool log of 1) timetep accepted, 2) 
      int itol;//parameter for scalar or vector tolerance computing
      std::vector<double> dt_v;
      std::vector<bool> accept_v;
      bool single_step_mode;
  };
}
#endif

