// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-10-10

#ifndef RL_H_IS_INCLUDED
#define RL_H_IS_INCLUDED

#include <boost/scoped_array.hpp>

#include "types.h"
#include "ODESolver.h"
#include "ODE.h"

namespace goss {

  // Forward declaration
  class LinearizedODE;

  // First order accurate Generalized Rush-Larsen ODE Solver
  class RL: public ODESolver
  {
  public:
    
    // Default Constructor
    RL();

    // Constructor
    RL(ODE* ode);

    // Copy constructor
    RL(const RL& solver);
    
    // Return a copy of itself
    ODESolver* copy() const { return new RL(*this); }

    // Destructor
    ~RL();

    // Attach ODE to solver
    virtual void attach(ODE* ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);
    
  private:

    // Storage for the Linearized ODE
    LinearizedODE* _lode;

    // Pointers to intermediate values used while stepping
    boost::scoped_array<double> a;
    boost::scoped_array<double> b;

    boost::scoped_array<uint> linear_terms;
  
  };

}
#endif
