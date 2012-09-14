#ifndef ODE_H_IS_INCLUDED
#define ODE_H_IS_INCLUDED

#include <vector>
#include <string>
#include <iostream>

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

    inline uint size() const { return _system_size; }
    virtual void eval(const double* state, double t, double* f_vals) = 0;
    virtual double eval(uint idx, const double* state, double t);

    virtual uint num_linear_terms() const {return 0;}
    virtual uint num_parameters() const {return 0;}
    virtual void linear_terms(int* indices) const {}
    virtual void linear_derivatives(const double* x, double t, double* y) const {}

    virtual void eval_intermediates(const double* x, double t, double* y)
    { /*TODO: Implement*/}
    
    virtual double eval_intermediate(uint i, const double* x, double t)
    { return 0.0; /*TODO: Implement*/}


    /* Method for extracting data from the ODE. Must be implemented in
     * subclasses. FIXME: Figure out a better system for this.
     */
    virtual void probe (double* y) { }

    virtual void set_parameters(const double* values) {}
    virtual void set_parameter(uint idx, double value) {}
    virtual double get_parameter(uint idx) const {}
    virtual void get_parameters(goss::DoubleVector*) const {}
    virtual std::string get_state_name(uint idx) const;
    virtual std::string get_parameter_name(uint idx) const;
    virtual void get_ic(goss::DoubleVector *res) const = 0;

  protected: 
    
    virtual void set_state_descriptions() {};
    virtual void set_parameter_names() {};
    
    const uint _system_size;
    const uint _parameter_size;
    
    std::vector<std::string> _state_descr;
    std::vector<std::string> _parameter_descr;
    
  };
}

#endif
