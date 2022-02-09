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

#ifndef ODESolver_h_IS_INCLUDED
#define ODESolver_h_IS_INCLUDED

#include <memory>
#include <iostream>

#include "types.h"
#include "ODE.h"
#include "Parameters.h"
#include "log.h"

//#include <fenv.h>

namespace goss
{

  class ODESolver {

  public:

    // Default parameters
    static Parameters default_parameters()
    {
      Parameters params("ode_solver");
      params.add("ldt", -1.0);
      return params;
    }

    // Default Constructor
    ODESolver () : parameters(), _ode(static_cast<ODE*>(0))
    {
      parameters = default_parameters();
    }

    // Copy constructor (Uses default copy constructor of ODE)
    ODESolver (const ODESolver& solver) : parameters(solver.parameters),
					  _ode(static_cast<ODE*>(0))
    {
      if (solver._ode)
	_ode = solver._ode->copy();
    }

    // Destructor
    virtual ~ODESolver ()
    {
      // Do nothing
    }

    // Return a copy of itself
    virtual std::shared_ptr<ODESolver> copy() const = 0;

    // Attach ODE and reset solver
    virtual void attach(std::shared_ptr<ODE> ode)
    {
      //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
      _ode = ode; reset();}

    // Reset solver
    virtual void reset() { /* Do nothing */ }

    // Step solver an interval of time forward
    virtual void forward(double* y, double t, double interval) = 0;

    // The size of the ODE
    inline uint num_states() const { return _ode ? _ode->num_states() : 0; }

    // Return the ODE (const version)
    inline const std::shared_ptr<ODE> get_ode() const { return _ode; }

    // Return the ODE
    inline std::shared_ptr<ODE> get_ode() { return _ode; }

    // Return the internal time step
    inline virtual double get_internal_time_step() const { return -1.; }

    // Set the internal time step
    inline virtual void set_internal_time_step(double) {}

    // Return true if the Solver is adaptive
    virtual bool is_adaptive() const { return false; }

    // Solver parameters
    Parameters parameters;

  protected:

    // Access to scratch space in ODE
    inline std::vector<double>& _f1() const
    {assert(_ode); return _ode->_f1;}

    inline std::vector<double>& _f2() const
    {assert(_ode); return _ode->_f2;}

    // Shared pointer to ode
    std::shared_ptr<ODE> _ode;

  };

}
#endif
