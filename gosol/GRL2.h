// Copyright (C) 2006 Skavhaug.
// All rights reserved.
//
// First added:  2007-07-09
// Last changed: 2008-11-10

#ifndef GRL2_H_IS_INCLUDED
#define GRL2_H_IS_INCLUDED

#include "ODESolver.h"
#include <math.h>


/* GRL2: Second order accurate Generalized Rush-Larsen ODE Solver
 *
 * A two step second order accurate ode solver.
 */

namespace gosol {

  class GRL2: public ODESolver
  {
      public:

        /// Constructor
        GRL2();
        GRL2(gosol::ODE* ode);
        virtual void attach(gosol::ODE* ode);
        void forward(double* y, double t, double dt);

        /// Destructor
        ~GRL2();

      private:
        double* y0; 
        double* a;
        double* b;
        int* linear_terms;
        const double delta;
        int n;
  };

}
#endif
