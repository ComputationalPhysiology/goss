#include <cmath>
#include <stdio.h>

#include "ThetaSolver.h"
#include "constants.h"
#include "log.h"

using namespace goss;

//-----------------------------------------------------------------------------
ThetaSolver::ThetaSolver() : ImplicitODESolver(), _z1(0), _ft1(0), _justrefined(false)
{
}
//-----------------------------------------------------------------------------
ThetaSolver::ThetaSolver(std::shared_ptr<ODE> ode)
    : ImplicitODESolver(), _z1(0), _ft1(0), _justrefined(false)
{
    attach(ode);
}
//-----------------------------------------------------------------------------
ThetaSolver::ThetaSolver(const ThetaSolver &solver)
    : ImplicitODESolver(solver), _z1(solver.num_states()), _ft1(solver.num_states()),
      _justrefined(solver._justrefined)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void ThetaSolver::attach(std::shared_ptr<ODE> ode)
{
    // Attach ode using bases
    ImplicitODESolver::attach(ode);

    // Init memory
    _z1.resize(num_states());
    _ft1.resize(num_states());
}
//-----------------------------------------------------------------------------
void ThetaSolver::reset()
{

    _justrefined = false;

    _jac_comp = 0;
    _stages = 1;

    ImplicitODESolver::reset();
}
//-----------------------------------------------------------------------------
void ThetaSolver::forward(double *y, double t, double dt)
{

    assert(_ode);

    uint i;
    bool step_ok;

    const double t_end = t + dt;
    const double ldt_0 = _ldt;

    double ldt = ldt_0 > 0 ? ldt_0 : dt;
    int num_refinements = 0;


    // A way to check if we are at t_end.
    const double eps = GOSS_EPS * 1000;

    for (i = 0; i < _ode->num_states(); ++i)
        _prev[i] = 0.0;

    while (true) {

        // Use 0.0 z1:
        for (i = 0; i < _ode->num_states(); ++i)
            _z1[i] = 0.0;

        // Explicit eval
        if (std::abs(theta - 1.0) > GOSS_EPS) {
            _ode->eval(y, t, &_ft1[0]);
            for (i = 0; i < _ode->num_states(); ++i)
                _prev[i] = (1 - theta) * _ft1[i];
        }

        // Check if we should re compute the jacobian
        if (num_refinements <= num_refinements_without_always_recomputing_jacobian)
            always_recompute_jacobian = true;

        // Solve for increment
        step_ok = newton_solve(&_z1[0], &_prev[0], y, t + theta * ldt, ldt, theta,
                               always_recompute_jacobian);

        // Newton step OK
        if (step_ok) {

            // Add increment
            for (i = 0; i < _ode->num_states(); ++i)
                y[i] += _z1[i];

            t += ldt;

            if (std::fabs(t - t_end) < eps)
                break;

            // If the solver has refined, we do not allow it to double its
            // timestep for anoter step
            if (!_justrefined) {

                // double time step
                const double tmp = 2.0 * ldt;
                if (ldt_0 > 0. && tmp >= ldt_0) {
                    ldt = ldt_0;
                } else {
                    ldt = tmp;
                    log(DBG, "Changing dt    | t : %g, from %g to %g", t, tmp / 2, ldt);
                }

            } else {
                _justrefined = false;
            }

            // If we are passed t_end
            if ((t + ldt + GOSS_EPS) > t_end) {
                ldt = t_end - t;
                log(DBG, "Changing ldt   | t : %g, to adapt for dt end: %g", t, ldt);
            }

        } else {
            ldt /= 2.0;
            if (ldt < min_dt) {
                goss_error("ThetaSolver.cpp", "Forward ThetaSolver",
                           "Newtons solver failed to converge as dt become smaller "
                           "than \"min_dt\" %e",
                           min_dt);
            }
            log(DBG, "Reducing dt    | t : %g, new: %g", t, ldt);
            _justrefined = true;
            num_refinements += 1;
        }
    }
#ifdef DEBUG
    // Lower level than DEBUG!
    log(5, "ThetaSolver done with comp_jac = %d and rejected = %d at t=%1.2e\n", _jac_comp,
        _rejects, t);
#endif
}
//-----------------------------------------------------------------------------
