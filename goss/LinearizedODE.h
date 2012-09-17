#ifndef LINEARIZED_ODE_H_IS_INCLUDED
#define LINEARIZED_ODE_H_IS_INCLUDED

#include "types.h"
#include "ODE.h"

namespace goss {

  // Class which provides an interface for ODE Solvers which need 
  // linearized terms
  class LinearizedODE : public virtual ODE
  {
  public:
    
    LinearizedODE(uint system_size) : ODE(system_size)
    { 
      // Do nothing
    } 

    virtual ~LinearizedODE() 
    {
      // Do nothing
    }

    // Populate indices with information about the linear terms
    virtual void linear_terms(uint* indices) const = 0;

    // Evaluate the linear derivatives
    virtual void linear_derivatives(const double* x, double t, double* y) const = 0;
    
  };
}

#endif
