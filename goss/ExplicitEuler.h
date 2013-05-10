// Copyright (C) 2006-2012 Ola Skavhaug
//
// This file is part of GOSS.
//
// GOSS is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GOSS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GOSS. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Johan Hake 2012

#ifndef ExplicitEuler_h_IS_INCLUDED
#define ExplicitEuler_h_IS_INCLUDED

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <cstdlib>

#include "ODESolver.h"

namespace goss 
{

  // An Explicit Euler solver
  class ExplicitEuler : public ODESolver
  {
    
    public:
    
    // Default constructor
    ExplicitEuler();

    // Constructor
    ExplicitEuler(double _ldt);
    
    // Constructor
    ExplicitEuler(boost::shared_ptr<ODE> ode, double _ldt=-1.0);

    // Copy constructor
    ExplicitEuler(const ExplicitEuler& solver);

    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const 
    { return boost::make_shared<ExplicitEuler>(*this); }

    // Attach ODE to solver
    virtual void attach(boost::shared_ptr<ODE> ode);

    // Step solver an interval in time forward
    void forward(double* y, double t, double interval);

    // Destructor
    ~ExplicitEuler ();
    
  protected:

    // State derivative, allocated in attach(ode)
    std::vector<double> _dFdt; 
    
  };
  
}
#endif
