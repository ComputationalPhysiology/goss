#ifndef AdaptiveImplicitSolver_h_IS_INCLUDED
#define AdaptiveImplicitSolver_h_IS_INCLUDED

#include "ImplicitODESolver.h"
#include <vector>
#include <stdio.h>

namespace goss 
{

  // Base class for adaptive implicit ODE solvers
  class AdaptiveImplicitSolver : public ImplicitODESolver
  {

  public:

    // Default constructor
    AdaptiveImplicitSolver();

    // Constructor
    AdaptiveImplicitSolver (double ldt);

    // Copy constructor
    AdaptiveImplicitSolver(const AdaptiveImplicitSolver& solver);

    // Destructor
    virtual ~AdaptiveImplicitSolver ();

    // Initialize data
    virtual void reset();
    
    // Step solver an interval of time forward
    virtual void forward(double* y, double t, double interval) = 0;

    // Return the current time
    double get_current_time();

    // Return the current time step
    double get_current_time_step();

    // Return number of accepted solutions
    long get_num_accepted(){return num_accepted;}

    // Return number of rejected solutions
    long get_num_rejected(){return num_rejected;}

    // FIXME: Where is this used!?
    // Store timestep and accepted timestep
    void log_data(double dt, bool accepted);
    
    // Return a vector of collected timesteps
    void dt_vector(DoubleVector *res);

    // Return a record of accepted 
    void accepted_vector(DoubleVector *res);

    // Set single step mode
    void set_single_step_mode(bool mode) {single_step_mode = mode;}

    // Set tolerance
    void set_tol(double atol, double rtol=1.0e-8) {_atol = atol; _rtol = rtol;}
    
    // Set iord
    void set_iord(int iord){_iord = iord;}

    // FIXME: Should this be protected?
    // Compute an initial time step guess
    double dtinit(double t, double* y0, double* y1, double* f0, double* f1, double iord);

    // FIXME: Should this be protected?
    // Compute new timestep 
    void new_time_step(double* y, double* yn, double* e, double t_end);

  protected: 
    
    // Log of 1) the numer of steps, 2) the number of rejected steps
    long num_accepted, num_rejected;

    // a bool log of 1) timetep accepted, 2) 
    bool step_accepted, reached_tend;

    std::vector<double> dt_v;
    std::vector<bool> accept_v;

    bool single_step_mode;
    
    double _t, _dt_prev;

    // local time step and tolerence.
    double _atol, _rtol, _iord, facmin, facmax, facmaxb, stabfac; 

    // Added stability to reduce the number of Jacobian computations
    double stabdown, stabup; 
    double err_old, dt_old;

    // Parameter for scalar or vector tolerance computing
    int _itol; 

  };

}
#endif

