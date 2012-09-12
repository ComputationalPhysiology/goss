// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2007-07-09

#ifndef RL_H_IS_INCLUDED
#define RL_H_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>


/// Documentation of class RL

namespace gossol {

  class RL: public ODESolver
  {
      public:

        /// Constructor
        RL();
        RL(gossol::ODE* ode);
        virtual void attach(gossol::ODE* ode);
        void forward(double* y, double t, double dt);

        /// Destructor
        ~RL();

      private:
        double* a;
        double* b;
        int* linear_terms;
        const double delta;
        int n;
  };

}
#endif
