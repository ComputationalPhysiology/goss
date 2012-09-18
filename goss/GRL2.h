// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-18

#ifndef GRL2_H_IS_INCLUDED
#define GRL2_H_IS_INCLUDED

#include "types.h"
#include "ODESolver.h"
#include "ODE.h"

namespace goss {

  // Forward declaration
  class LinearizedODE;

  // Second order accurate Generalized Rush-Larsen ODE Solver
  class GRL2 : public ODESolver
  {
  public:

    // Default Constructor
    GRL2();

    // Constructor
    GRL2(ODE* ode);
    
    // Destructor
    ~GRL2();

    // Attach ODE to solver
    virtual void attach(ODE* ode);

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
