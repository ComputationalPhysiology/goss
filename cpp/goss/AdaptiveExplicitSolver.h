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

#ifndef AdaptiveExplicitSolver_h_IS_INCLUDED
#define AdaptiveExplicitSolver_h_IS_INCLUDED

#include "ODESolver.h"
#include <vector>

namespace goss
{

// Base class for Adaptive and Expicit solvers
class AdaptiveExplicitSolver : public ODESolver
{

  public:
    // Default Constructor
    AdaptiveExplicitSolver();

    // Copy constructor
    AdaptiveExplicitSolver(const AdaptiveExplicitSolver &solver);

    // Destructor
    virtual ~AdaptiveExplicitSolver(){};

    // Initialize Solver
    virtual void reset();

    // Return the current time
    double get_current_time();

    // Return the current time step
    double get_current_time_step();

    // Return number of accepted solutions
    long get_num_accepted()
    {
        return num_accepted;
    }

    // Return number of rejected solutions
    long get_num_rejected()
    {
        return num_rejected;
    }

    // Set single step mode
    void set_single_step_mode(bool mode)
    {
        single_step_mode = mode;
    }

    // Set tolerance
    void set_tol(double atol, double rtol = 1.0e-8)
    {
        _atol = atol;
        _rtol = rtol;
    }

    // Get the absolute tolerance
    double get_atol()
    {
        return _atol;
    }

    // Get the relative tolerance
    double get_rtol()
    {
        return _rtol;
    }

    // Set iord
    void set_iord(int iord)
    {
        _iord = iord;
    }

    // Get iord
    int get_iord()
    {
        return _iord;
    }

    // FIXME: Should this be protected?
    // Compute an initial time step guess
    double dtinit(double t, double *y0, double *y1, double *f0, double *f1, double iord);

    // FIXME: Should this be protected?
    // Compute new timestep
    void new_time_step(double *y, double *yn, double *e, double t_end);

    // Return true if the Solver is adaptive
    bool is_adaptive() const
    {
        return true;
    }


  protected:
    // Log of 1) the numer of steps, 2) the number of rejected steps
    ulong num_accepted, num_rejected;
    double _t, _ldt, _dt, _dt_prev;

    // Local time step and tolerence.
    double _atol, _rtol, _iord, facmin, facmax, facmaxb, stabfac;

    // Bool log of 1) timestep accepted, 2)
    bool step_accepted, reached_tend;

    // Parameter for scalar or vector tolerance computing
    int _itol;
    std::vector<double> dt_v;
    std::vector<bool> accept_v;
    bool single_step_mode;
};
} // namespace goss
#endif
