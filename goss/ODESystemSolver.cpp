#include <iostream>
#include <stdexcept>

#ifdef HAS_OPENMP
#include <omp.h>
#endif

#include "ODESystemSolver.h"

using namespace goss;

//-----------------------------------------------------------------------------
ODESystemSolver::ODESystemSolver(uint num_nodes, ODESolver* solver, 
				 ParameterizedODE* ode) :
  _num_nodes(num_nodes), _num_threads(0), _solver(solver), 
  _threaded_solvers(0), _ode(ode), _states(num_nodes*ode->num_states()),
  _field_parameters(num_nodes*ode->num_field_parameters()),
  _ldt_vec(solver->is_adaptive() ? num_nodes : 0, solver->get_internal_time_step()),
  _is_adaptive(solver->is_adaptive()), 
  _has_field_parameters(ode->num_field_parameters() > 0)
{ 
  

  // Attach ODE to solver
  _solver->attach(ode);

  // Reset values for the field parameters and the states
  reset_default();
  
}
//-----------------------------------------------------------------------------
ODESystemSolver::~ODESystemSolver()
{ 
  for (uint i = 0; i < _threaded_solvers.size(); i++)
    delete _threaded_solvers[i];
}
//-----------------------------------------------------------------------------
void ODESystemSolver::forward(double t, double interval)
{

  // Local solver pointer
  ODESolver* solver = _solver.get();

  // Iterate over all nodes using threaded or non-threaded loop
  if(_num_threads > 0)
  {
#pragma omp parallel for schedule(guided, 20) 
   for (uint node = 0; node < _num_nodes; node++)
      _forward_node(_threaded_solvers[omp_get_thread_num()], node, t, interval);
  }
  else
  {
    for (uint node = 0; node < _num_nodes; node++)
      _forward_node(solver, node, t, interval);
  }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::get_field_states(double* values, bool tangled_storage) const
{
  
  // Iterate over all nodes using threaded or non-threaded loop
  if(_num_threads > 0)
  {
#pragma omp parallel for schedule(guided, 20) 
   for (uint node = 0; node < _num_nodes; node++)
     _get_field_states_node(node, values, tangled_storage);
  }
  else
  {
    for (uint node = 0; node < _num_nodes; node++)
     _get_field_states_node(node, values, tangled_storage);
  }
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_field_states(const double* values, bool tangled_storage)
{
  
  // Iterate over all nodes using threaded or non-threaded loop
  if(_num_threads > 0)
  {
#pragma omp parallel for schedule(guided, 20) 
   for (uint node = 0; node < _num_nodes; node++)
     _set_field_states_node(node, values, tangled_storage);
  }
  else
  {
    for (uint node = 0; node < _num_nodes; node++)
     _set_field_states_node(node, values, tangled_storage);
  }
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_field_parameters(const double* values, bool tangled_storage)
{

  
  // Iterate over all nodes using threaded or non-threaded loop
  if(_num_threads > 0)
  {
#pragma omp parallel for schedule(guided, 20) 
   for (uint node = 0; node < _num_nodes; node++)
     _set_field_parameters_node(node, values, tangled_storage);
  }
  else
  {
    for (uint node = 0; node < _num_nodes; node++)
     _set_field_parameters_node(node, values, tangled_storage);
  }
  
}
//-----------------------------------------------------------------------------
void ODESystemSolver::reset_default()
{
  std::string param;

  // Get default ic of ODE
  DoubleVector default_ic;
  _ode->get_ic(&default_ic);

  // Get default field parameters
  std::vector<double> default_field_params(_ode->num_field_parameters(), 0.0);
  for (uint i = 0; i < _ode->num_field_parameters(); i++)
  {
    param = _ode->get_field_parameter_names()[i];
    default_field_params[i] = _ode->get_parameter(param);
  }

  // Iterate over nodes and set default values
  if(_num_threads > 0)
  {
#pragma omp parallel for schedule(guided, 20) 
    for (uint node = 0; node < _num_nodes; node++)
      _reset_default_node(node, default_ic, default_field_params);
  }
  else
  {
    for (uint node = 0; node < _num_nodes; node++)
      _reset_default_node(node, default_ic, default_field_params);
  }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_num_threads(uint num_threads)
{
#ifdef HAS_OPENMP
  // Store num threads
  _num_threads = num_threads;

  // Set num threads in Open MP library
  omp_set_num_threads(num_threads);

  // Delete and resize
  for (uint i = 0; i < _threaded_solvers.size(); i++)
    delete _threaded_solvers[i];
  _threaded_solvers.resize(num_threads);

  // Re-create threaded solvers
  for (uint i = 0; i<num_threads; i++)
    _threaded_solvers[i] = _solver->copy();
  
#else
  
  // Keeps the compiler happy
  _num_threads = num_threads*0;
#endif
}
//-----------------------------------------------------------------------------
