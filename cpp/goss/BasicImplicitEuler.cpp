#include <cassert>
#include <cmath>
#include <cstdio>

#include "BasicImplicitEuler.h"
#include "constants.h"
#include "log.h"


using namespace goss;

//-----------------------------------------------------------------------------
BasicImplicitEuler::BasicImplicitEuler() : ImplicitODESolver()
{
}
//-----------------------------------------------------------------------------
BasicImplicitEuler::BasicImplicitEuler(std::shared_ptr<ODE> ode) : ImplicitODESolver()
{
    attach(ode);
}
//-----------------------------------------------------------------------------
BasicImplicitEuler::BasicImplicitEuler(const BasicImplicitEuler &solver) : ImplicitODESolver(solver)
{

    _stages = 1;
}
//-----------------------------------------------------------------------------
BasicImplicitEuler::~BasicImplicitEuler()
{

    _stages = 1;
}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::attach(std::shared_ptr<ODE> ode)
{

    // Attach ode using bases
    ImplicitODESolver::attach(ode);
}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::reset()
{

    ImplicitODESolver::reset();
}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::compute_factorized_jacobian(double *y, double t, double dt)
{

    // Let ODE compute the jacobian
    _ode->compute_jacobian(y, t, _jac.data());

    // Build Euler discretization of jacobian
    mult(-dt, _jac.data());
    add_mass_matrix(_jac.data());

    // Factorize the jacobian
    _ode->lu_factorize(_jac.data());
    _jac_comp += 1;
}
//-----------------------------------------------------------------------------
void BasicImplicitEuler::forward(double *y, double t, double dt)
{
    // FIXME: This does not work on Mac?
    // feenableexcept(FE_INVALID | FE_OVERFLOW);
    assert(_ode);

    // Calculate number of steps and size of timestep based on _ldt
    const double ldt_0 = _ldt;
    const ulong nsteps = ldt_0 > 0 ? std::ceil(dt / ldt_0 - 1.0E-12) : 1;
    const double ldt = dt / nsteps;

    const double rtol = relative_tolerance;

    // Local time
    double lt = t;

    // Calculate number of steps and size of timestep based on _ldt

    /* USE THIS WHEN WE HAVE ONE NEWTON ITERATION INPLACE
  const ulong nsteps = _ldt > 0 ? std::ceil(dt/_ldt - 1.0E-12) : 1;
  const double dt = dt/nsteps;

  for (ulong step = 0; step < nsteps; ++step)
  {
    // Evaluate rhs
    _ode->eval(y, lt, &_dFdt[0]);

    // Update states
    for (uint i = 0; i < num_states(); ++i)
      y[i] += dt*_dFdt[i];

    // Increase time
    lt += dt;
  }

  for (i = 0; i < num_states(); ++i)
    _prev[i] = 0.0;
  */

    for (ulong step = 0; step < nsteps; ++step) {

        double relative_residual = 1.0;
        double initial_residual = 1.0;
        double residual;

        bool newton_converged = false;
        _newton_iterations = 0;

        // Copy previous solution
        for (uint i = 0; i < num_states(); ++i)
            _prev[i] = y[i];

        // Start iterations
        while (!newton_converged && _newton_iterations < max_iterations) {

            // Evaluate ODE using computed solution
            _ode->eval(y, lt + ldt, _f1.data());

            // Build rhs for linear solve
            // b = eval(y) - (y-y0)/dt
            for (uint i = 0; i < num_states(); ++i)
                _b[i] = (y[i] - _prev[i]) * _ode->differential_states()[i] - dt * _f1[i];

            // Compute residual
            residual = norm(_b.data());

            if (_newton_iterations == 0)
                initial_residual = residual;

            // Relative residual
            relative_residual = residual / initial_residual;

            if (relative_residual < rtol) {
                newton_converged = true;
                break;
            }

            // Compute Jacobian
            compute_factorized_jacobian(y, lt + ldt, ldt);

            // Linear solve on factorized jacobian
            _ode->forward_backward_subst(_jac.data(), _b.data(), _dz.data());

            // Compute initial residual
            if (_newton_iterations == 0)
                initial_residual = residual;

            for (uint i = 0; i < num_states(); ++i)
                y[i] -= _dz[i];

            // Update number of iterations
            ++_newton_iterations;

            // Output iteration number and residual
            log(DBG,
                "BasicImplicitEuler newton iteration %d: r (abs) = %.3e "
                " r (rel) = %.3e (tol = %.3e)",
                _newton_iterations, residual, relative_residual, rtol);
        }

        if (!newton_converged)
            error("Newton solver did not converge. Maximal newton iterations exceded.");
        else
            // Output iteration number and residual
            log(DBG,
                "BasicImplicitEuler newton iteration %d: r (abs) = %.3e "
                " r (rel) = %.3e (tol = %.3e)",
                _newton_iterations, residual, relative_residual, rtol);

        // Increase local time
        lt += dt;
    }

#ifdef DEBUG
    // Lower level than DEBUG!
    log(5, "BasicImplicitEuler done with comp_jac = %d at t=%1.2e\n", _jac_comp, t);
#endif
}
//-----------------------------------------------------------------------------
