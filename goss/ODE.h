#ifndef ODE_H_IS_INCLUDED
#define ODE_H_IS_INCLUDED

#include <vector>
#include <string>

#include "types.h"
#include "DoubleVector.h"

namespace goss {

  class ODE 
  {
  public:
    ODE(uint system_size, uint parameter_size) : 
      _system_size(system_size), _parameter_size(parameter_size), 
      _state_descr(system_size, ""), _parameter_descr(parameter_size, "")
    { 
      // Do nothing
    } 

    virtual ~ODE() 
    {
      // Do nothing
    }

    // Return the size of the ODE
    inline uint size() const { return _system_size; }

    // Evaluate rhs of the ODE
    virtual void eval(const double* state, double t, double* f_vals) = 0;

    // Evaluate component idx of the rhs of the ODE
    virtual double eval(uint idx, const double* state, double t);

    // The number of linear terms in the rhs
    virtual uint num_linear_terms() const {return 0;}

    // The number of parameters
    virtual uint num_parameters() const {return 0;}

    // Populate indices with information about the linear terms
    virtual void linear_terms(int* indices) const {}

    // Evaluate the linear derivatives
    virtual void linear_derivatives(const double* x, double t, double* y) const {}

    // Evaluate the intermediates
    virtual void eval_intermediates(const double* x, double t, double* y)
    { /*TODO: Implement*/}
    
    // Evaluate componentwise intermediates
    virtual double eval_intermediate(uint i, const double* x, double t)
    { return 0.0; /*TODO: Implement*/}

    // Set all parameters
    virtual void set_parameters(const double* values) {}

    // Set single parameter
    virtual void set_parameter(uint idx, double value) {}

    // Get parameter idx
    virtual double get_parameter(uint idx) const {}

    // Get all parameters
    virtual void get_parameters(goss::DoubleVector*) const {}

    // Get state name idx
    virtual std::string get_state_name(uint idx) const;

    // Get parameter name idx
    virtual std::string get_parameter_name(uint idx) const;

    // Get default initial conditions
    virtual void get_ic(goss::DoubleVector *res) const = 0;

  protected: 
    
    // Set state descriptions
    virtual void set_state_descriptions() {};
    
    // Set parameter names
    virtual void set_parameter_names() {};
    
    // ODE size
    const uint _system_size;

    // Number of parameters
    const uint _parameter_size;
    
    // State descriptions
    std::vector<std::string> _state_descr;
    
    // Parameter descriptions
    std::vector<std::string> _parameter_descr;
    
  };
}

#endif
