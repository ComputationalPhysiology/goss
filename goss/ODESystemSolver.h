#ifndef ODESYSTEMSOLVER_h_IS_INCLUDED
#define ODESYSTEMSOLVER_h_IS_INCLUDED

#include <vector>
#include <boost/scoped_ptr.hpp>

#include "types.h"
#include "DoubleVector.h"

namespace goss 
{

  class ParameterizedODE;
  class ODESolver;

  // class to handle the update systems of the same ODE
  class ODESystemSolver
  {
    
  public:
    
    // Constructor
    ODESystemSolver(uint nodes, ODESolver* solver, ParameterizedODE* ode);
    
    // Step all nodes an interval of time forward
    void forward(double t, double interval);

    // Return field state values
    void get_field_states(double* values, bool tangled_storage=true) const;

    // Set field state values
    void set_field_states(const double* values, bool tangled_storage=true);
    
    // Set field parameter values
    void set_field_parameters(const double* values, bool tangled_storage);
    
    // Use initial condition and initial values of field parameters to 
    // reset System variables
    void reset_default();

    // Set the number of threads
    inline void set_num_threads(uint num_threads);
    
    // Get the number of threads
    inline uint get_num_threads();
    
  private:

    // Help functions used in normal and OpenMP runs
    // ---------------------------------------------

    // Reset default values for a given node
    inline void _reset_default_node(uint node, const DoubleVector& default_ic,
				    const std::vector<double>& default_field_params);

    // Forward states for a given node
    inline void _forward_node(ODESolver* solver, uint node, double t, double interval);

    // Set field states for a given node
    inline void _set_field_states_node(uint node, const double* values, 
				       bool tangled_storage);

    // Set field states for a given node
    inline void _get_field_states_node(uint node, double* values, 
				       bool tangled_storage) const;

    // Set field parameter values
    inline void _set_field_parameters_node(uint node, const double* values, 
					   bool tangled_storage);

    // Number of nodes
    const uint _num_nodes;

    // The ODE solver
    boost::scoped_ptr<ODESolver> _solver;

    // Local pointer to ODE (No ownership. Solver owes ODE)
    ParameterizedODE* _ode;
    
    // Solution array with the solution for all nodes
    std::vector<double> _states;
    
    // Field parameters
    std::vector<double> _field_parameters;
    
    // Storage of local time step for each node
    std::vector<double> _ldt_vec;

    // Flag for having adaptive ODE solver
    const bool _is_adaptive;

    // Flag for having field parameters
    const bool _has_field_parameters;

  };
  
}
#endif
