#ifndef ImplicitODESolver_h_IS_INCLUDED
#define ImplicitODESolver_h_IS_INCLUDED

#include "ODESolver.h"
#include <iostream>
#include <math.h>
#include <cstdlib>

namespace gosol 
{

  class ImplicitODESolver :public ODESolver {
    public:
      void setNewtonTol(double newton_tol_=1e-5){newton_tol=newton_tol_;}

      ~ImplicitODESolver () {
        if (jac != NULL) {
          for (int i=0;i< jac_size; ++i)
            free(jac[i]);
        }
        free(jac);
        free(f1);
        free(f2);
        free(yz);
        free(b); 
        free(dz);
        free(prev);
      }

      ImplicitODESolver();
      ImplicitODESolver(double ldt);

    protected: 

      //variables used in the jacobian evaluation
      double **jac;
      double *f1, *f2, *yz;
      double dt_old; //h; 

      double *b, *dz; // right hand side and solution of the linear system
      double *prev; //previous stages, used by DIRK methods
      double newton_tol;
      double eta;// Variable used in the estimation of the error of the newton iteration for the first iteration (Very important for linear problems!)
      double kappa;// The safety factor for the stopping criterion of the newton iteration
      int jac_size;

      int stages;
      int newtonits, maxits, rejects, jac_comp;
      double  min_dt;
      bool recompute_jacobian;

      void computeJacobian(double t, double* y);
      void mult(double fact, double** matrix);
      void addIdentity(double** matrix);

      void factLU(double** mat);
      void forwBackLU(const double * const * mat, double* b, double* x);

      void init(); 


      //this function is designed for SDIRK and Backward Euler:
      virtual bool NewtonSolve(double* k, double* prev, double* y0, double t, double dt, double alpha);
      virtual double norm(double* vec);
  };

}
#endif
