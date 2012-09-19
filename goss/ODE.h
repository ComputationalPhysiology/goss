#ifndef ODE_H_IS_INCLUDED
#define ODE_H_IS_INCLUDED

#include "types.h"
#include "DoubleVector.h"

namespace goss {

  class ODE 
  {
  public:

    // Base class for an ODE
    ODE(uint num_states_) : 
      _num_states(num_states_)
    { 
      // Do nothing
    } 

    virtual ~ODE() 
    {
      // Do nothing
    }

    // Return the size of the ODE
    inline uint num_states() const { return _num_states; }

    // Evaluate rhs of the ODE
    virtual void eval(const double* state, double t, double* f_vals) = 0;

    // Evaluate component idx of the rhs of the ODE
    virtual double eval(uint idx, const double* state, double t);

    // Get default initial conditions
    virtual void get_ic(goss::DoubleVector *res) const = 0;

    // Return a copy of the ODE
    virtual ODE* copy() const = 0;

  protected: 
    
    // ODE size
    const uint _num_states;

  };
}

#endif
