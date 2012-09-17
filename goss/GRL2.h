// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-17

#ifndef GRL2_H_IS_INCLUDED
#define GRL2_H_IS_INCLUDED

#include "types.h"
#include "ODESolver.h"
#include "LinearizedODE.h"

namespace goss {

  // Second order accurate Generalized Rush-Larsen ODE Solver
  class GRL2 : public ODESolver
  {
  public:

    // Default Constructor
    GRL2();

    // Constructor
    GRL2(goss::LinearizedODE* ode);
    
    // Destructor
    ~GRL2();

    // Attach ODE to solver
    virtual void attach(LinearizedODE* ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);
    
  private:

    // Storage for the Linearized ODE
    LinearizedODE* _lode;

    // Pointers to intermediate values used while stepping
    double* y0; 
    double* a;
    double* b;

    uint* linear_terms;
    const double delta;

    // Number of bytes which will be copied each time step
    uint nbytes;

  };
  
}
#endif
