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

#ifndef ImplicitEuler_h_IS_INCLUDED
#define ImplicitEuler_h_IS_INCLUDED

#include <memory>
#include <vector>

#include "ImplicitODESolver.h"

namespace goss
{

// Implicit Euler
class ImplicitEuler : public ImplicitODESolver
{
  public:
    // Default parameters
    int num_refinements_without_always_recomputing_jacobian = 2;
    double min_dt = 0.0001;

    // Default constructor
    ImplicitEuler();

    // Constructor
    ImplicitEuler(std::shared_ptr<ODE> ode);

    // Copy constructor
    ImplicitEuler(const ImplicitEuler &solver);

    // Return a copy of itself
    std::shared_ptr<ODESolver> copy() const
    {
        return std::make_shared<ImplicitEuler>(*this);
    }

    // Destructor
    ~ImplicitEuler();

    // Attach ODE
    void attach(std::shared_ptr<ODE> ode);

    // Reset ODE
    void reset();

    // Step solver an interval of time forward
    void forward(double *y, double t, double interval);

  protected:
    std::vector<double> _z1;
    bool _justrefined;
};

} // namespace goss
#endif
