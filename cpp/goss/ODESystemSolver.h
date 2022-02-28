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

#ifndef ODESYSTEMSOLVER_h_IS_INCLUDED
#define ODESYSTEMSOLVER_h_IS_INCLUDED

#include <memory>
#include <vector>

#include "DoubleVector.h"
#include "ODESolver.h"
#include "ParameterizedODE.h"
#include "log.h"
#include "types.h"

bool has_openmp();

namespace goss
{

class ODESolver;

// class to handle the update systems of the same ODE
class ODESystemSolver
{

  public:
    // Constructor
    ODESystemSolver(uint nodes, std::shared_ptr<ODESolver> solver,
                    std::shared_ptr<ParameterizedODE> ode);

    // Destructor
    ~ODESystemSolver();

    // Step all nodes an interval of time forward
    void forward(double t, double interval);

    // Step all nodes an multiple time steps
    void solve(double *field_states, double *t, const ulong num_timesteps,
               bool tangled_storage = true);

    // Return components of system field state values
    void get_field_state_components(double *component_field_states, uint num_components,
                                    const uint *components, bool tangled_storage = true) const;

    // Set components of system field state values
    void set_field_state_components(const double *component_field_states, uint num_components,
                                    const uint *components, bool tangled_storage = true);

    // Return system field state values
    void get_field_states(double *system_field_states, bool tangled_storage = true) const;

    // Set system field state values
    void set_field_states(const double *system_field_states, bool tangled_storage = true);

    // Reutrn system field parameter values
    void get_field_parameters(double *system_field_params, bool tangled_storage = true) const;

    // Set system field parameter values
    void set_field_parameters(const double *system_field_params, bool tangled_storage = true);

    // Use initial condition and initial values of field parameters to
    // reset System variables
    void reset_default();

    // Set the number of threads
    void set_num_threads(uint num_threads);

    // Get a const pointer to the whole states data array
    const double *states() const
    {
        return _states.data();
    }

    // Get a pointer to the whole states data array
    double *states()
    {
        return _states.data();
    }

    // Get a const pointer to the states data of one node
    const double *states(uint node) const
    {
        return &_states[node * _ode->num_states()];
    }

    // Get a pointer to the states data of one node
    double *states(uint node)
    {
        return &_states[node * _ode->num_states()];
    }

    // Get the number of threads
    inline uint get_num_threads() const;

    // Get the ParameterizedODE
    inline std::shared_ptr<ParameterizedODE> ode()
    {
        return _ode;
    }

    // Get the ODESolver
    inline std::shared_ptr<ODESolver> solver()
    {
        return _solver;
    }

    // Get the number of nodes
    inline uint num_nodes()
    {
        return _num_nodes;
    }

  private:
    // Help functions used in normal and OpenMP runs
    // ---------------------------------------------

    // Reset default values for a given node
    inline void _reset_default_node(uint node, const DoubleVector &default_ic,
                                    const std::vector<double> &default_field_params);

    // Forward states for a given node
    inline void _forward_node(ODESolver &solver, uint node, double t, double interval);

    // Set field states for a given node
    inline void _set_field_states_node(uint node, const double *values, bool tangled_storage);

    // Get field states for a given node
    inline void _get_field_states_node(uint node, double *values, bool tangled_storage) const;


    // Set field state components for a given node
    inline void _set_field_states_node_comp(uint node, const double *values, uint num_components,
                                            const uint *components, bool tangled_storage);

    // Get field state components for a given node
    inline void _get_field_states_node_comp(uint node, double *values, uint num_components,
                                            const uint *components, bool tangled_storage) const;

    // Get field parameter values
    inline void _get_field_parameters_node(uint node, double *values, bool tangled_storage) const;

    // Set field parameter values
    inline void _set_field_parameters_node(uint node, const double *values, bool tangled_storage);

    // Number of nodes
    const uint _num_nodes;

    // Number of threads
    uint _num_threads;

    // The ODE solver
    std::shared_ptr<ODESolver> _solver;

    // Solvers used in OpenMP threaded runs
    std::vector<std::shared_ptr<ODESolver>> _threaded_solvers;

    // Local pointer to ODE (No ownership. Solver owes ODE)
    std::shared_ptr<ParameterizedODE> _ode;

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
void ODESystemSolver::_reset_default_node(uint node, const DoubleVector &default_ic,
                                          const std::vector<double> &default_field_params)
{

    // Set field parameters
    if (_has_field_parameters) {

        // Update field parameters
        for (uint i = 0; i < _ode->num_field_parameters(); i++) {
            _field_parameters[_ode->num_field_parameters() * node + i] = default_field_params[i];
        }
    }

    // Set states
    for (uint i = 0; i < _ode->num_states(); i++)
        _states[_ode->num_states() * node + i] = default_ic.data[i];
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_forward_node(ODESolver &solver_, uint node, double t, double interval)
{

    // Update field parameters if any
    if (_has_field_parameters) {

        // Get local ode. Nead to grab ode from solver so it also works in a
        // threader version
        ParameterizedODE &ode_ = dynamic_cast<ParameterizedODE &>(*solver_.get_ode());
        ode_.set_field_parameters(&_field_parameters[node * ode_.num_field_parameters()]);
    }

    // Set internal time step if adaptive
    if (_is_adaptive)
        solver_.set_internal_time_step(_ldt_vec[node]);

    // Forward solver_
    solver_.forward(&_states[node * _ode->num_states()], t, interval);

    // Get internal time step if adaptive
    if (_is_adaptive)
        _ldt_vec[node] = solver_.get_internal_time_step();
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_set_field_states_node(uint node, const double *values, bool tangled_storage)
{
    uint ind = 0;

    for (uint i = 0; i < _ode->num_field_states(); i++) {
        ind = tangled_storage ? node * _ode->num_field_states() + i : _num_nodes * i + node;
        _states[node * _ode->num_states() + _ode->get_field_state_indices()[i]] = values[ind];
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_get_field_states_node_comp(uint node, double *values, uint num_components,
                                                  const uint *components,
                                                  bool tangled_storage) const
{
    uint ind = 0;
    uint i = 0;

    for (uint ic = 0; ic < num_components; ic++) {
        i = components[ic];
        ind = tangled_storage ? node * _ode->num_field_states() + i : _num_nodes * i + node;

        values[ind] = _states[node * _ode->num_states() + _ode->get_field_state_indices()[i]];
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_set_field_states_node_comp(uint node, const double *values,
                                                  uint num_components, const uint *components,
                                                  bool tangled_storage)
{
    uint ind = 0;
    uint i = 0;

    for (uint ic = 0; ic < num_components; ic++) {
        i = components[ic];
        ind = tangled_storage ? node * _ode->num_field_states() + i : _num_nodes * i + node;
        _states[node * _ode->num_states() + _ode->get_field_state_indices()[i]] = values[ind];
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_get_field_states_node(uint node, double *values, bool tangled_storage) const
{
    uint ind = 0;

    for (uint i = 0; i < _ode->num_field_states(); i++) {
        ind = tangled_storage ? node * _ode->num_field_states() + i : _num_nodes * i + node;

        values[ind] = _states[node * _ode->num_states() + _ode->get_field_state_indices()[i]];
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_get_field_parameters_node(uint node, double *values,
                                                 bool tangled_storage) const
{
    uint ind = 0;

    for (uint i = 0; i < _ode->num_field_parameters(); i++) {
        ind = tangled_storage ? node * _ode->num_field_parameters() + i : _num_nodes * i + node;

        values[ind] = _field_parameters[node * _ode->num_field_parameters() + i];
    }
}
//-----------------------------------------------------------------------------
void ODESystemSolver::_set_field_parameters_node(uint node, const double *values,
                                                 bool tangled_storage)
{

    uint ind = 0;

    for (uint i = 0; i < _ode->num_field_parameters(); i++) {
        ind = tangled_storage ? node * _ode->num_field_parameters() + i : _num_nodes * i + node;
        _field_parameters[node * _ode->num_field_parameters() + i] = values[ind];
    }
}
//-----------------------------------------------------------------------------
} // namespace goss
#endif
