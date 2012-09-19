#ifndef ImplicitODESolver_h_IS_INCLUDED
#define ImplicitODESolver_h_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include "ODESolver.h"

namespace goss 
{

  // Base class of all Implicit Solvers
  class ImplicitODESolver : public ODESolver 
  {
  
  public:
  
    // Default constructor
    ImplicitODESolver();
  
    // Constructor
    ImplicitODESolver(double ldt);

    // Copy constructor
    ImplicitODESolver(const ImplicitODESolver& solver);

    // Destructor
    virtual ~ImplicitODESolver ();

    // Reset solver
    virtual void reset();

    // Attach ode
    virtual void attach(ODE* ode);

    // Step solver an interval of time forward
    virtual void forward(double* y, double t, double interval) = 0;

    // Set newton tolerance
    void set_newton_tol(double newton_tol){ _newton_tol = newton_tol; }

    // Return the number of recomputation of the jacobian
    int num_jac_comp(){ return jac_comp; }

  protected: 

    // Compute jacobian using numerical approximation
    void compute_jacobian(double t, double* y);

    // Scale a matrix
    void mult(double fact, boost::scoped_array<boost::scoped_array<double> >& matrix);

    // Add identity to matrix
    void add_identity(boost::scoped_array<boost::scoped_array<double> >& matrix);

    // LU Factorize matrix
    void lu_factorize(boost::scoped_array<boost::scoped_array<double> >& matrix);

    // Forward/Backward supstituion of factories matrix
    void forward_backward_subst(
            const boost::scoped_array<boost::scoped_array<double> >& mat, 
	    double* b, double* x);

    // This function is designed for SDIRK and Backward Euler:
    virtual bool newton_solve(double* k, double* prev, double* y0, double t, 
			      double dt, double alpha);

    // Compute the norm of a vector
    virtual double norm(double* vec);

    // Variables used in the jacobian evaluation
    boost::scoped_array<boost::scoped_array<double> > jac;
    boost::scoped_array<double> f1, f2, yz;

    // Right hand side and solution of the linear system
    boost::scoped_array<double> _b, dz;

    // Previous stages, used by DIRK methods
    boost::scoped_array<double> _prev;

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
    ulong num_tsteps;

    double min_dt;
    bool recompute_jacobian;

  };

}
#endif
