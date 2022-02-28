// Copyright (C) 2012 Johan Hake
//
// This file is part of GOSS.
//
// GOSS is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GOSS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GOSS. If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
bool has_openmp()
{
    return true;
}
#else
bool has_openmp()
{
    return false;
}
int omp_get_thread_num()
{
    return 0;
}
#endif

#include "ODESystemSolver.h"

using namespace goss;

//-----------------------------------------------------------------------------
ODESystemSolver::ODESystemSolver(uint num_nodes_, std::shared_ptr<ODESolver> solver_,
                                 std::shared_ptr<ParameterizedODE> ode_)
    : _num_nodes(num_nodes_), _num_threads(0), _solver(solver_), _threaded_solvers(0), _ode(ode_),
      _states(num_nodes_ * ode_->num_states()),
      _field_parameters(num_nodes_ * ode_->num_field_parameters()),
      _ldt_vec(solver_->is_adaptive() ? num_nodes_ : 0, solver_->get_internal_time_step()),
      _is_adaptive(solver_->is_adaptive()), _has_field_parameters(ode_->num_field_parameters() > 0)
{

    // Attach ODE to solver
    _solver->attach(ode_);

    // Reset values for the field parameters and the states
    reset_default();
}
//-----------------------------------------------------------------------------
ODESystemSolver::~ODESystemSolver()
{
    //  for (uint i = 0; i < _threaded_solvers.size(); i++)
    //    delete _threaded_solvers[i];
}
//-----------------------------------------------------------------------------
void ODESystemSolver::forward(double t, double interval)
{

    // Iterate over all nodes using threaded or non-threaded loop
    if (_num_threads > 0) {
#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _forward_node(*_threaded_solvers[omp_get_thread_num()], node, t, interval);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _forward_node(*_solver, node, t, interval);
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::solve(double *field_states, double *t, const ulong num_timesteps,
                            bool tangled_storage)
{

    double t_next, t_current;
    double dt;
    uint num_field_states = _ode->num_field_states();

    t_current = t[0];
    ulong i, j, it;
    for (it = 1; it < num_timesteps; it++) {
        t_next = t[it];
        dt = t_next - t_current;
        forward(t_current, dt);
        get_field_states(field_states + it * _num_nodes * num_field_states, tangled_storage);

        t_current = t_next;
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::get_field_state_components(double *component_field_states,
                                                 uint num_components, const uint *components,
                                                 bool tangled_storage) const
{
    // Iterate over all nodes using threaded or non-threaded loop
    if (_num_threads > 0) {
        //#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _get_field_states_node_comp(node, component_field_states, num_components, components,
                                        tangled_storage);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _get_field_states_node_comp(node, component_field_states, num_components, components,
                                        tangled_storage);
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_field_state_components(const double *component_field_states,
                                                 uint num_components, const uint *components,
                                                 bool tangled_storage)
{
    // Iterate over all nodes using threaded or non-threaded loop
    if (_num_threads > 0) {
        //#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _set_field_states_node_comp(node, component_field_states, num_components, components,
                                        tangled_storage);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _set_field_states_node_comp(node, component_field_states, num_components, components,
                                        tangled_storage);
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::get_field_states(double *field_state_values, bool tangled_storage) const
{

    // Iterate over all nodes using threaded or non-threaded loop
    if (_num_threads > 0) {
        //#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _get_field_states_node(node, field_state_values, tangled_storage);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _get_field_states_node(node, field_state_values, tangled_storage);
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_field_states(const double *field_state_values, bool tangled_storage)
{

    // Iterate over all nodes using threaded or non-threaded loop
    if (_num_threads > 0) {
        //#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _set_field_states_node(node, field_state_values, tangled_storage);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _set_field_states_node(node, field_state_values, tangled_storage);
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::get_field_parameters(double *field_param_values, bool tangled_storage) const
{

    // Iterate over all nodes using threaded or non-threaded loop
    if (_num_threads > 0) {
        //#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _get_field_parameters_node(node, field_param_values, tangled_storage);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _get_field_parameters_node(node, field_param_values, tangled_storage);
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_field_parameters(const double *field_param_values, bool tangled_storage)
{


    // Iterate over all nodes using threaded or non-threaded loop
    if (_num_threads > 0) {
        //#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _set_field_parameters_node(node, field_param_values, tangled_storage);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _set_field_parameters_node(node, field_param_values, tangled_storage);
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
    for (uint i = 0; i < _ode->num_field_parameters(); i++) {
        param = _ode->get_field_parameter_names()[i];
        default_field_params[i] = _ode->get_parameter(param);
    }

    // Iterate over nodes and set default values
    if (_num_threads > 0) {
#pragma omp parallel for schedule(guided, 20)
        for (uint node = 0; node < _num_nodes; node++)
            _reset_default_node(node, default_ic, default_field_params);
    } else {
        for (uint node = 0; node < _num_nodes; node++)
            _reset_default_node(node, default_ic, default_field_params);
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::set_num_threads(uint num_threads)
{
#ifdef _OPENMP

    char *GOMP_CPU_AFFINITY;
    GOMP_CPU_AFFINITY = getenv("GOMP_CPU_AFFINITY");
    if (GOMP_CPU_AFFINITY != NULL)
        printf("GOMP_CPU_AFFINITY: %s\n", GOMP_CPU_AFFINITY);

    // Store num threads
    _num_threads = num_threads;

    // Set num threads in Open MP library
    omp_set_num_threads(num_threads);

    // Delete and resize
    //for (uint i = 0; i < _threaded_solvers.size(); i++)
    //  delete _threaded_solvers[i];
    _threaded_solvers.resize(num_threads);

    // Re-create threaded solvers
    for (uint i = 0; i < num_threads; i++)
        _threaded_solvers[i] = _solver->copy();

#else

    // Keeps the compiler happy
    _num_threads = num_threads * 0;

#endif
}
//-----------------------------------------------------------------------------
