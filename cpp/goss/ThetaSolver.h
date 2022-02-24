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

#include <memory>
#include <vector>

#include "ImplicitODESolver.h"

namespace goss
{

class ThetaSolver : public ImplicitODESolver
{
  public:
    // Default parameters
    int num_refinements_without_always_recomputing_jacobian = 2;
    double min_dt = 0.0001;
    double theta = 0.5;

    // Default constructor
    ThetaSolver();

    // Constructor
    ThetaSolver(std::shared_ptr<ODE> ode);

    // Copy constructor
    ThetaSolver(const ThetaSolver &solver);

    // Destructor
    ~ThetaSolver(){};

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const
    {
        return std::make_shared<ThetaSolver>(*this);
    }

    // Attach ODE
    void attach(std::shared_ptr<ODE> ode);

    // Reset solver
    void reset();

    // Step solver an interval of time forward
    void forward(double *y, double t, double dt);

  protected:
    std::vector<double> _z1, _ft1;
    bool _justrefined;
};

} // namespace goss
#endif
