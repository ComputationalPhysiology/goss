#ifndef PARAMETERIZED_ODE_H_IS_INCLUDED
#define PARAMETERIZED_ODE_H_IS_INCLUDED

#include <vector>
#include <string>

#include "types.h"
#include "ODE.h"
#include "DoubleVector.h"

namespace goss {

  // Class which provides a more verbose interface for users to controll an ODE
  class ParameterizedODE : public virtual ODE
  {
  public:

    // Constructor
    ParameterizedODE(uint system_size, uint num_parameters_, uint num_field_states_, 
		     uint num_field_parameters_, uint num_intermediates_) : 
      ODE(system_size), 
      _state_names(system_size, ""), 
      _field_state_names(num_field_states_, ""), 
      _parameter_names(num_parameters_, ""), 
      _field_parameter_names(num_field_parameters_, ""), 
      _field_state_indices(num_field_states_, 0), 
      _intermediate_names(num_intermediates_, "")
    { 
      // Do nothing
    } 

    virtual ~ParameterizedODE() 
    {
      // Do nothing
    }

    // Return the number of field states
    inline uint num_field_states() const { return _field_state_names.size(); }

    // The number of parameters
    inline uint num_parameters() const { return _parameter_names.size(); }

    // The number of field parameters
    inline uint num_field_parameters() const { return _field_parameter_names.size(); }

    // The number of field parameters
    inline uint num_intermediates() const { return _intermediate_names.size(); }

    // Evaluate the intermediates
    virtual void eval_intermediates(const double* x, double t, double* y) const = 0;
    
    // Evaluate componentwise intermediates
    virtual double eval_intermediate(uint i, const double* x, double t) const = 0;

    // Set all field parameters
    virtual void set_field_parameters(const double* values) = 0;

    // Set a parameter
    virtual void set_parameter(std::string name, double value) = 0;

    // Get a parameter value
    virtual double get_parameter(std::string name) const = 0;

    // Get all state names
    const std::vector<std::string>& get_state_names() const
    { return _state_names; }

    // Get field state names
    const std::vector<std::string>& get_field_state_names() const
    { return _field_state_names; }

    // Get all parameter names
    const std::vector<std::string>& get_parameter_names() const
    { return _parameter_names; }

    // Get field parameter names
    const std::vector<std::string>& get_field_parameter_names() const
    { return _field_parameter_names; }

    // Get field parameter names
    const std::vector<uint>& get_field_state_indices() const
    { return _field_state_indices; }

    // Get intermediate names
    const std::vector<std::string>& get_intermediate_names() const
    { return _intermediate_names; }

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

    // Vector with intermediate names
    std::vector<std::string> _intermediate_names;
    
  };
}

#endif
