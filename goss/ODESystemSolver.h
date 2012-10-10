#ifndef ODESYSTEMSOLVER_h_IS_INCLUDED
#define ODESYSTEMSOLVER_h_IS_INCLUDED

#include <vector>
#include <boost/scoped_ptr.hpp>

#include "types.h"
#include "DoubleVector.h"
#include "ParameterizedODE.h"
#include "ODESolver.h"

namespace goss 
{

  class ODESolver;

  // class to handle the update systems of the same ODE
  class ODESystemSolver
  {
    
  public:
    
    // Constructor
    ODESystemSolver(uint nodes, ODESolver* solver, ParameterizedODE* ode);
    
    // Destructor
    ~ODESystemSolver();
    
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
    void set_num_threads(uint num_threads);
    
    // Get the number of threads
    inline uint get_num_threads() const; 
    
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

    // Number of threads
    uint _num_threads;

    // The ODE solver
    boost::scoped_ptr<ODESolver> _solver;

    // Solvers used in OpenMP threaded runs
    std::vector<ODESolver* > _threaded_solvers;

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
  //---------------------------------------------------------------------------
  uint ODESystemSolver::get_num_threads() const
  {
    return _num_threads;
  }
  //---------------------------------------------------------------------------
  void ODESystemSolver::_reset_default_node(uint node, const DoubleVector& default_ic,
					  const std::vector<double>& default_field_params)
  {
  
    // Set field parameters
    if (_has_field_parameters)
    {
      
      // Update field parameters
      for (uint i = 0; i < _ode->num_field_parameters(); i++)
      {
	_field_parameters[_ode->num_field_parameters()*node + i] =	\
	  default_field_params[i];
      }
    }
    
    // Set states
    for (uint i = 0; i < _ode->num_states(); i++)
      _states[_ode->num_states()*node + i] = default_ic.data[i];
    
  }
  //-----------------------------------------------------------------------------
  void ODESystemSolver::_forward_node(ODESolver* solver, uint node, double t, double interval)
  {
    
  // Update field parameters if any
    if (_has_field_parameters)
    {
      
      // Get local ode. Nead to grab ode from solver so it also works in a 
      // threader version
    ParameterizedODE* ode = dynamic_cast<ParameterizedODE*>(solver->get_ode());
    ode->set_field_parameters(&_field_parameters[node*ode->num_field_parameters()]);
    }
    
    // Set internal time step if adaptive
    if (_is_adaptive)
      solver->set_internal_time_step(_ldt_vec[node]);
    
    // Forward solver
    solver->forward(&_states[node*_ode->num_states()], t, interval);
    
    // Get internal time step if adaptive
    if (_is_adaptive)
      _ldt_vec[node] = solver->get_internal_time_step();
    
  }
  //-----------------------------------------------------------------------------
  void ODESystemSolver::_set_field_states_node(uint node, const double* values, 
					       bool tangled_storage)
  {
    uint ind = 0;
    
    for (uint i = 0; i < _ode->num_field_states(); i++)
    {
      ind = tangled_storage ? node*_ode->num_field_states() + i :	\
	_num_nodes*i + node;
      _states[node*_ode->num_states() + _ode->get_field_state_indices()[i]] = values[ind];
    }
    
  }
//-----------------------------------------------------------------------------
  void ODESystemSolver::_get_field_states_node(uint node, double* values, 
					     bool tangled_storage) const
  {
    uint ind = 0;
    
    for (uint i = 0; i < _ode->num_field_states(); i++)
    {
      ind = tangled_storage ? node*_ode->num_field_states() + i :	\
	_num_nodes*i + node;
      
      values[ind] = _states[node*_ode->num_states() +		\
			    _ode->get_field_state_indices()[i]];
    }
    
  }
//-----------------------------------------------------------------------------
  void ODESystemSolver::_set_field_parameters_node(uint node, const double* values, 
						   bool tangled_storage)
  {
    
    uint ind = 0;
    
    for (uint i = 0; i < _ode->num_field_parameters(); i++)
    {
      ind = tangled_storage ? node*_ode->num_field_parameters() + i :	\
	_num_nodes*i + node;
      _field_parameters[node*_ode->num_field_parameters() +		\
			_ode->get_field_state_indices()[i]] = values[ind];
    }
    
  }
//-----------------------------------------------------------------------------
}
#endif

