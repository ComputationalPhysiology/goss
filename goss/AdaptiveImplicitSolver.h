#ifndef AdaptiveImplicitSolver_h_IS_INCLUDED
#define AdaptiveImplicitSolver_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <math.h>
#include <vector>
#include <stdio.h>

namespace goss 
{

  class AdaptiveImplicitSolver : public ImplicitODESolver
  {

    public:
      AdaptiveImplicitSolver (goss::ODE* ode_);
      AdaptiveImplicitSolver();
      ~AdaptiveImplicitSolver ();

      void init();
      /*virtual void attach(goss::ODE* ode_);*/
      virtual void forward(double* y, double t, double interval) = 0;

      double getCurrentTime();
      double getCurrentTimeStep();
      long getNumAccepted(){return num_accepted;}
      long getNumRejected(){return num_rejected;}
      void logData(double dt, bool accepted);
      void dtVector(goss::DoubleVector *res);
      void acceptedVector(goss::DoubleVector *res);
      int numJacComp(){return jac_comp;}

      void setSingleStepMode(bool mode){single_step_mode=mode;}
      void setTol(double atol_, double rtol_=1.0e-8);
      void setIord(int iord_);
      double dtinit(double t, double* y0, double* y1, double* f0, double* f1, double iord);
      void newTimeStep(double* y, double* yn, double* e, double t_end);

    protected: 
      long num_accepted, num_rejected;// a log of 1) the numer of steps, 2)the number of rejected steps
      bool step_accepted, reached_tend;// a bool log of 1) timetep acceptec, 2) 
      //#ifdef DEBUG
      std::vector<double> dt_v;
      std::vector<bool> accept_v;
      //#endif
      bool single_step_mode;

      double t, dt_prev;
      double atol, rtol, iord, facmin, facmax, facmaxb, stabfac; // local time step and tolerence.
      double stabdown, stabup; // Added stability to reduce the number if Jacobian computations
      double err_old, dt_old;
      int itol; // Parameter for scalar or vector tolerance computing
  };

}
#endif

