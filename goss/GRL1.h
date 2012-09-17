// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-09-17

#ifndef GRL1_H_IS_INCLUDED
#define GRL1_H_IS_INCLUDED

#include "types.h"
#include "ODESolver.h"
#include "LinearizedODE.h"

namespace goss {

  // First order accurate Generalized Rush-Larsen ODE Solver
  class GRL1: public ODESolver
  {
  public:
    
    // Default Constructor
    GRL1();

    // Constructor
    GRL1(goss::LinearizedODE* ode);
    
    // Destructor
    ~GRL1();

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
    const double delta;

  };
}
#endif
