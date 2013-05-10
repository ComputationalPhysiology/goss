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

#ifndef ThetaSolver_h_IS_INCLUDED
#define ThetaSolver_h_IS_INCLUDED

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "ImplicitODESolver.h"

namespace goss 
{

  class ThetaSolver: public ImplicitODESolver
  {
  public:
  
    // Default constructor
    ThetaSolver() {};
    
    // Constructor
    ThetaSolver(double ldt); 
    
    // Constructor
    ThetaSolver(boost::shared_ptr<ODE> ode, double ldt=-1.0);
    
    // Destructor
    ~ThetaSolver(){};

    // Return a copy of itself
    boost::shared_ptr<ODESolver> copy() const 
    { return boost::make_shared<ThetaSolver>(*this); }

    // Attach ODE
    void attach(boost::shared_ptr<ODE> ode);

    // Reset solver
    void reset();

    // Step solver an interval of time forward
    void forward(double* y, double t, double dt);

    // Theta parameter
    double theta;

    std::vector<int> newton_iter1;
    std::vector<int> newton_accepted1;
    std::vector<double> dt_v;
    
    protected:
    
    std::vector<double> z1, ft1;
    bool justrefined;
    
  };
  
}
#endif
