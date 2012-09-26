#include <iostream>
#include <stdexcept>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

#include "ODESystemSolver.h"
#include "ParameterizedODE.h"
#include "ODESolver.h"

using namespace goss;

//-----------------------------------------------------------------------------
ODESystemSolver::ODESystemSolver(uint num_nodes, ODESolver* solver, 
				 ParameterizedODE* ode) :
  _num_nodes(num_nodes), _solver(solver), _ode(ode), 
  _states(num_nodes*ode->num_states()),
  _field_parameters(num_nodes*ode->num_field_parameters()),
  _ldt_vec(solver->is_adaptive() ? num_nodes : 0, solver->get_internal_time_step()),
  _is_adaptive(solver->is_adaptive()), 
  _has_field_parameters(ode->num_field_parameters() > 0)
{ 

  // Reset values for the field parameters and the states
  reset_default();
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::forward(double t, double interval)
{

  // Local solver pointer
  ODESolver* solver = _solver.get();
  const uint num_threads = get_num_threads();

  // Iterate over all nodes
  // FIXME: Move if tests outside update loop?
  #pragma omp parallel for if(num_threads > 0) schedule(guided, 20) firstprivate(solver)
  for (uint node = 0; node < _num_nodes; node++)
    _forward_node(solver, node, t, interval);

}
//-----------------------------------------------------------------------------
void ODESystemSolver::get_field_states(double* values, bool tangled_storage) const
{
  
  // FIXME: Add logics for OpenMP
  for (uint node = 0; node < _num_nodes; node++)
    _get_field_states_node(node, values, tangled_storage);
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_field_states(const double* values, bool tangled_storage)
{
  
  // FIXME: Add logics for OpenMP
  for (uint node = 0; node < _num_nodes; node++)
    _set_field_states_node(node, values, tangled_storage);
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_field_parameters(const double* values, bool tangled_storage)
{

  // FIXME: Add logics for OpenMP
  for (uint node = 0; node < _num_nodes; node++)
    _set_field_parameters_node(node, values, tangled_storage);
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::reset_default()
{
  std::string param;

  // Get default ic of ODE
  DoubleVector default_ic;
  _ode->get_ic(&default_ic);

  // Get default field parameters
  std::vector<double> default_field_params(_ode->num_field_states(), 0.0);
  for (uint i = 0; i < _ode->num_field_parameters(); i++)
  {
    param = _ode->get_field_parameter_names()[i];
    default_field_params[i] = _ode->get_parameter(param);
  }

  // Iterate over nodes and set default values
  for (uint node = 0; node < _num_nodes; node++)
    _reset_default_node(node, default_ic, default_field_params);
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_num_threads(uint num_threads)
{
#ifdef HAS_OPENMP
  omp_set_num_threads(num_threads);
#else
  numthreads;
#endif
}
//-----------------------------------------------------------------------------
uint ODESystemSolver::get_num_threads()
{
#ifdef HAS_OPENMP
  return omp_get_num_threads();
#else
  return 0;
#endif
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_reset_default_node(uint node, const DoubleVector& default_ic,
					  const std::vector<double>& default_field_params)
{
  
  // Set field parameters
  if (_has_field_parameters)
  {
    
    // Update field parameters
    for (uint i = 0; i < _ode->num_field_parameters(); i++)
    {
      _field_parameters[_ode->num_field_parameters()*node + i] = \
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
    ind = tangled_storage ? node*_ode->num_field_states() + i : \
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
    ind = tangled_storage ? node*_ode->num_field_states() + i : \
      _num_nodes*i + node;
    values[ind] = _states[node*_ode->num_states() + \
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
    ind = tangled_storage ? node*_ode->num_field_parameters() + i : \
      _num_nodes*i + node;
    _field_parameters[node*_ode->num_field_parameters() + \
		      _ode->get_field_state_indices()[i]] = values[ind];
  }
  
}
//-----------------------------------------------------------------------------

