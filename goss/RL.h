// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-17

#ifndef RL_H_IS_INCLUDED
#define RL_H_IS_INCLUDED

#include "types.h"
#include "ODESolver.h"
#include "LinearizedODE.h"

namespace goss {

  // First order accurate Generalized Rush-Larsen ODE Solver
  class RL: public ODESolver
  {
  public:
    
    // Default Constructor
    RL();

    // Constructor
    RL(goss::LinearizedODE* ode);
    
    // Destructor
    ~RL();

    // Attach ODE to solver
    virtual void attach(goss::LinearizedODE* ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);
    
  private:

    // Storage for the Linearized ODE
    LinearizedODE* _lode;

    // Pointers to intermediate values used while stepping
    double* a;
    double* b;

    uint* linear_terms;
  
  };

}
#endif
