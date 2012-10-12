// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2012-10-11

#ifndef GRL1_H_IS_INCLUDED
#define GRL1_H_IS_INCLUDED

#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "types.h"
#include "ODESolver.h"
#include "ODE.h"

namespace goss {

  // Forward declaration
  class LinearizedODE;

  // First order accurate Generalized Rush-Larsen ODE Solver
  class GRL1: public ODESolver
  {
  public:
    
    // Default Constructor
    GRL1();

    // Constructor
    GRL1(boost::shared_ptr<ODE> ode);
    
    // Copy constructor
    GRL1(const GRL1& solver);
    
    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const { return boost::make_shared<GRL1>(*this); }

    // Destructor
    ~GRL1();

    // Attach ODE to solver
    virtual void attach(boost::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);
    
  private:

    // Store a downcasted pointer to the Linearized ODE
    LinearizedODE* _lode;

    // Pointers to intermediate values used while stepping
    boost::scoped_array<double> a;
    boost::scoped_array<double> b;

    boost::scoped_array<uint> linear_terms;
    const double delta;

  };
}
#endif
