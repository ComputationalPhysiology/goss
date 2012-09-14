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
    ParameterizedODE(uint system_size, uint parameter_size) : 
      ODE(system_size), _parameter_size(parameter_size)
    { 
      // Do nothing
    } 

    virtual ~ParameterizedODE() 
    {
      // Do nothing
    }

    // The number of parameters
    inline uint num_parameters() const {return _parameter_size;}

    // Evaluate the intermediates
    virtual void eval_intermediates(const double* x, double t, double* y) const = 0;
    
    // Evaluate componentwise intermediates
    virtual double eval_intermediate(uint i, const double* x, double t) const = 0;

    // Set all parameters
    virtual void set_parameters(const double* values) = 0;

    // Set single parameter
    virtual void set_parameter(uint idx, double value) = 0;

    // Get parameter idx
    virtual double get_parameter(uint idx) const = 0;

    // Get all parameters
    virtual void get_parameters(goss::DoubleVector*) const = 0;

    // Get state name idx
    virtual std::string get_state_name(uint idx) const;

    // Get parameter name idx
    virtual std::string get_parameter_name(uint idx) const;

  protected: 
    
    // FIXME: Force user to fill these vectors in constructor?
    // Set state descriptions
    virtual void set_state_descriptions() = 0;
    
    // Set parameter names
    virtual void set_parameter_names() = 0;
    
    // Number of parameters
    const uint _parameter_size;
    
    // State descriptions
    std::vector<std::string> _state_descr;
    
    // Parameter descriptions
    std::vector<std::string> _parameter_descr;
    
  };
}

#endif
