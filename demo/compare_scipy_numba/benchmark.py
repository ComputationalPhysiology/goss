import timeit
from pathlib import Path

import goss
import matplotlib.pyplot as plt
import numba
import numpy as np
from scipy.integrate import solve_ivp


def main(rhs, odefile, t, args, run_plot=True, run_timings=True):  # noqa: C901
    numba_rhs = numba.jit(rhs, nopython=True)

    ode = goss.ParameterizedODE(odefile)
    name = Path(odefile).stem
    u0 = ode.get_ic()

    tspan = [t[0], t[-1]]

    goss_ExplicitEuler = goss.solvers.ExplicitEuler(ode)
    goss_ExplicitEuler.internal_time_step = 0.001

    goss_RL1 = goss.solvers.RL1(ode)
    goss_RL1.internal_time_step = 0.001

    goss_GRL1 = goss.solvers.GRL1(ode)
    goss_GRL1.internal_time_step = 0.001

    def _solve_scipy(f, method):
        return solve_ivp(
            f,
            tspan,
            y0=u0,
            args=args,
            method=method,
            rtol=1e-8,
            atol=1e-8,
            t_eval=t,
        ).y

    def solve_scipy_RK45():
        return _solve_scipy(rhs, method="RK45")

    def solve_scipy_numba_RK45():
        return _solve_scipy(numba_rhs, method="RK45")

    def solve_scipy_RK23():
        return _solve_scipy(rhs, method="RK23")

    def solve_scipy_numba_RK23():
        return _solve_scipy(numba_rhs, method="RK23")

    def solve_goss_ExplicitEuler():
        return goss_ExplicitEuler.solve(t)

    def solve_goss_RL1():
        return goss_RL1.solve(t)

    def solve_goss_GRL1():
        return goss_GRL1.solve(t)

    def time_methods():
        number = 1
        repeat = 2

        time_outfile = Path(f"timings_{name}_{number}_{repeat}.npy")

        if not time_outfile.is_file():
            print("Time methods")
            timings = {}

            for label, solver in [
                ("goss (ExplicitEuler)", solve_goss_ExplicitEuler),
                ("goss (RL1)", solve_goss_RL1),
                ("goss (GRL1)", solve_goss_GRL1),
                ("scipy (RK45)", solve_scipy_RK45),
                ("scipy+numba (RK45)", solve_scipy_numba_RK45),
                ("scipy (RK23)", solve_scipy_RK23),
                ("scipy+numba (RK23)", solve_scipy_numba_RK23),
            ]:
                print(f"Run benchmark for {label}")
                times = timeit.Timer(solver).repeat(repeat=repeat, number=number)
                timings[label] = min(times) / number
                np.save(time_outfile, timings)

        timings = np.load(time_outfile, allow_pickle=True).item()
        x = np.arange(len(timings))
        fig, ax = plt.subplots()
        ax.bar(x, timings.values(), align="center", log=True)
        ax.set_xticks(x)
        ax.set_xticklabels(timings.keys(), rotation=30)
        ax.set_ylabel("Time [seconds]")
        ax.set_xlabel("Solver")
        ax.grid()

        fig.savefig("timings_lorentz.png")
        plt.show()

    def plot_methods():

        solutions_outfile = Path(f"solutions_{name}.npy")

        if not solutions_outfile.is_file():
            data = {}
            for label, solver in [
                ("goss (ExplicitEuler)", solve_goss_ExplicitEuler),
                ("goss (RL1)", solve_goss_RL1),
                ("goss (GRL1)", solve_goss_GRL1),
                ("scipy (RK45)", solve_scipy_RK45),
                ("scipy+numba (RK45)", solve_scipy_numba_RK45),
                ("scipy (RK23)", solve_scipy_RK23),
                ("scipy+numba (RK23)", solve_scipy_numba_RK23),
            ]:
                print(f"Solve {label}")
                y = solver()
                if y.shape[1] == len(t):
                    y = y.T

                data[label] = y
            np.save(solutions_outfile, data)

        data = np.load(solutions_outfile, allow_pickle=True).item()

        fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
        lines = []
        labels = []
        for label, y in data.items():
            (l,) = ax[0].plot(t, y[:, 0])
            ax[1].plot(t, y[:, 1])
            ax[2].plot(t, y[:, 2])
            lines.append(l)
            labels.append(label)

        ax[0].set_ylabel("$x$")
        ax[1].set_ylabel("$z$")
        ax[2].set_ylabel("$z$")
        ax[2].set_xlabel("$t$")

        lgd = fig.legend(lines, labels, title="Solver", loc="center right")
        fig.subplots_adjust(right=0.75)
        fig.savefig(
            f"solutions_{name}.png",
            bbox_extra_artists=(lgd,),
            bbox_inches="tight",
            dpi=300,
        )
        plt.show()

    if run_plot:
        plot_methods()

    if run_timings:
        time_methods()
