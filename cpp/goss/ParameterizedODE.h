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

#ifndef PARAMETERIZED_ODE_H_IS_INCLUDED
#define PARAMETERIZED_ODE_H_IS_INCLUDED

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "DoubleVector.h"
#include "ODE.h"
#include "types.h"

namespace goss
{

// Class which provides a more verbose interface for users to controll an ODE
class ParameterizedODE : public ODE
{
  public:
    // Constructor
    ParameterizedODE(uint num_states_, uint num_parameters_, uint num_field_states_,
                     uint num_field_parameters_, uint num_monitored_)
        : ODE(num_states_), _state_names(num_states_, ""),
          _field_state_names(num_field_states_, ""), _parameter_names(num_parameters_, ""),
          _field_parameter_names(num_field_parameters_, ""),
          _field_state_indices(num_field_states_, 0), _monitored_names(num_monitored_, ""),
          _param_to_value()
    {
        // Do nothing
    }

    // Copy constructor
    ParameterizedODE(const ParameterizedODE &ode)
        : ODE(ode), _state_names(ode._state_names), _field_state_names(ode._field_state_names),
          _parameter_names(ode._parameter_names),
          _field_parameter_names(ode._field_parameter_names),
          _field_state_indices(ode._field_state_indices), _monitored_names(ode._monitored_names),
          _param_to_value(ode._param_to_value)
    {
        // Do nothing
    }

    virtual ~ParameterizedODE()
    {
        // Do nothing
    }

    // Return the number of field states
    inline uint num_field_states() const
    {
        return _field_state_names.size();
    }

    // The number of parameters
    inline uint num_parameters() const
    {
        return _parameter_names.size();
    }

    // The number of field parameters
    inline uint num_field_parameters() const
    {
        return _field_parameter_names.size();
    }

    // The number of field parameters
    inline uint num_monitored() const
    {
        return _monitored_names.size();
    }

    // Evaluate the monitored
    virtual void eval_monitored(const double *states, double t, double *monitored) const = 0;

    // Evaluate monitored values for many time steps
    virtual void monitored_values(const double *states, const double *t, double *monitored,
                                  double *m, const ulong num_timesteps) const
    {
        ulong i, j;
        double ti;
        for (i = 0; i < num_timesteps; i++) {
            ti = t[i];
            eval_monitored(states + i * num_states(), ti, m);
            for (j = 0; j < num_monitored(); j++) {
                monitored[i * num_monitored() + j] = m[j];
            }
        }
    }

    // Set all field parameters
    virtual void set_field_parameters(const double *field_params) = 0;

    // Set a parameter
    void set_parameter(std::string name, double value)
    {
        *_find_parameter(name) = value;
    }

    // Get a parameter value
    double get_parameter(std::string name) const
    {
        return *_find_parameter(name);
    }

    // Get all state names
    const std::vector<std::string> &get_state_names() const
    {
        return _state_names;
    }

    // Get field state names
    const std::vector<std::string> &get_field_state_names() const
    {
        return _field_state_names;
    }

    // Get all parameter names
    const std::vector<std::string> &get_parameter_names() const
    {
        return _parameter_names;
    }

    // Get field parameter names
    const std::vector<std::string> &get_field_parameter_names() const
    {
        return _field_parameter_names;
    }

    // Get field state indices
    const std::vector<uint> &get_field_state_indices() const
    {
        return _field_state_indices;
    }

    // Get monitored names
    const std::vector<std::string> &get_monitored_names() const
    {
        return _monitored_names;
    }

  protected:
    // These vectors are initialized in this class, but need to be
    // filled in derived classes.

    // Vector with state names
    std::vector<std::string> _state_names;

    // Vector with field state names
    std::vector<std::string> _field_state_names;

    // Vector with all parameter names
    std::vector<std::string> _parameter_names;

    // Vector with field parameter names
    std::vector<std::string> _field_parameter_names;

    // Vector with field state indices
    std::vector<uint> _field_state_indices;

    // Vector with monitored names
    std::vector<std::string> _monitored_names;

    // A map between parameter name and its value
    std::map<std::string, double *> _param_to_value;

  private:
    // Find a parameter value
    inline double *_find_parameter(std::string name) const
    {
        // Search parameter map
        std::map<std::string, double *>::const_iterator it = _param_to_value.find(name);
        if (it == _param_to_value.end())
            throw std::runtime_error("\"" + name + "\" is not a parameter in  "
				 "this ODE.");

        // Return value pointer
        return it->second;
    }
};
} // namespace goss

#endif
