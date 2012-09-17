#ifndef ImplicitODESolver_h_IS_INCLUDED
#define ImplicitODESolver_h_IS_INCLUDED

#include "ODESolver.h"

namespace goss 
{

  // Base class of all Implicit Solvers
  class ImplicitODESolver : public ODESolver 
  {
  
  public:
  
    // Constructor
    ImplicitODESolver();
  
    // Constructor
    ImplicitODESolver(double ldt);

    // Destructor
    virtual ~ImplicitODESolver ();

    // Set newton tolerance
    void set_newton_tol(double newton_tol){_newton_tol = newton_tol;}

    // Return the number of recomputation of the jacobian
    int num_jac_comp(){return jac_comp;}

  protected: 

    // Variables used in the jacobian evaluation
    double **jac;
    double *f1, *f2, *yz;
    double dt_old;

    // Right hand side and solution of the linear system
    double *_b, *dz;

    // Previous stages, used by DIRK methods
    double *_prev;

    // Newton tolerance
    double _newton_tol;

    // Variable used in the estimation of the error of the newton 
    // iteration for the first iteration (Important for linear problems!)
    double eta;

    // The safety factor for the stopping criterion of the newton iteration
    double kappa;
    uint jac_size;

    uint stages;
    int newtonits, maxits, rejects, jac_comp;
    double  min_dt;
    bool recompute_jacobian;

    // Compute jacobian using numerical approximation
    void compute_jacobian(double t, double* y);

    // Scale a matrix
    void mult(double fact, double** matrix);

    // Add identity to matrix
    void add_identity(double** matrix);

    // LU Factorize matrix
    void lu_factorize(double** mat);

    // Forward/Backward supstituion of factories matrix
    void forward_backward_subst(const double* const* mat, double* b, double* x);

    // Init solver
    void init(); 

    // This function is designed for SDIRK and Backward Euler:
    virtual bool newton_solve(double* k, double* prev, double* y0, double t, 
			      double dt, double alpha);
    virtual double norm(double* vec);
  };

}
#endif
