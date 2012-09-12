// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2007-07-09

#ifndef GRL1_H_IS_INCLUDED
#define GRL1_H_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>


/// Documentation of class GRL1

namespace gosol {

  class GRL1: public ODESolver
  {
      public:

        /// Constructor
        GRL1();
        GRL1(gosol::ODE* ode);
        virtual void attach(gosol::ODE* ode);
        void forward(double* y, double t, double dt);

        /// Destructor
        ~GRL1();

      private:
        double* a;
        double* b;
        int* linear_terms;
        const double delta;
        int n;
  };

}
#endif
