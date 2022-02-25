import timeit
from pathlib import Path

import goss
import matplotlib.pyplot as plt
import numba
import numpy as np
from scipy.integrate import solve_ivp


def default_internal_time_steps():
    return {
        "goss (ExplicitEuler)": -1,
        "goss (RK2)": -1,
        "goss (RK4)": -1,
        "goss (RL1)": -1,
        "goss (RL2)": -1,
        "goss (GRL1)": -1,
        "goss (GRL2)": -1,
        "goss (ImplicitEuler)": -1,
        "goss (ThetaSolver)": -1,
        "goss (RKF32)": -1,
        "goss (ESDIRK23a)": -1,
    }


def main(  # noqa: C901
    rhs,
    odefile,
    t,
    args,
    selected_methods=None,
    selected_states=None,
    run_plot=True,
    run_timings=True,
    recompute=False,
    internal_time_step=None,
    number=10,
    repeat=10,
):
    numba_rhs = numba.jit(rhs, nopython=True)
    internal_time_steps = default_internal_time_steps()
    if internal_time_step is not None:
        internal_time_steps.update(internal_time_step)

    ode = goss.ParameterizedODE(odefile)
    name = Path(odefile).stem
    u0 = ode.get_ic()
    if selected_states is None:
        selected_states = np.arange(len(u0))

    tspan = [t[0], t[-1]]

    goss_ExplicitEuler = goss.solvers.ExplicitEuler(ode)
    goss_ExplicitEuler.internal_time_step = internal_time_steps["goss (ExplicitEuler)"]

    goss_RK2 = goss.solvers.RK2(ode)
    goss_RK2.internal_time_step = internal_time_steps["goss (RK2)"]

    goss_RK4 = goss.solvers.RK4(ode)
    goss_RK4.internal_time_step = internal_time_steps["goss (RK4)"]

    goss_RL1 = goss.solvers.RL1(ode)
    goss_RL1.internal_time_step = internal_time_steps["goss (RL1)"]

    goss_RL2 = goss.solvers.RL2(ode)
    goss_RL2.internal_time_step = internal_time_steps["goss (RL2)"]

    goss_GRL1 = goss.solvers.GRL1(ode)
    goss_GRL1.internal_time_step = internal_time_steps["goss (GRL1)"]

    goss_GRL2 = goss.solvers.GRL2(ode)
    goss_GRL2.internal_time_step = internal_time_steps["goss (GRL2)"]

    goss_ImplicitEuler = goss.solvers.ImplicitEuler(ode)
    goss_ImplicitEuler.internal_time_step = internal_time_steps["goss (ImplicitEuler)"]

    goss_ThetaSolver = goss.solvers.ThetaSolver(ode)
    goss_ThetaSolver.internal_time_step = internal_time_steps["goss (ThetaSolver)"]

    goss_RKF32 = goss.solvers.RKF32(ode)
    goss_RKF32.internal_time_step = internal_time_steps["goss (RKF32)"]

    goss_ESDIRK23a = goss.solvers.ESDIRK23a(ode)
    goss_ESDIRK23a.internal_time_step = internal_time_steps["goss (ESDIRK23a)"]

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

    def solve_scipy_DOP853():
        return _solve_scipy(rhs, method="DOP853")

    def solve_scipy_numba_DOP853():
        return _solve_scipy(numba_rhs, method="DOP853")

    def solve_scipy_Radau():
        return _solve_scipy(rhs, method="Radau")

    def solve_scipy_numba_Radau():
        return _solve_scipy(numba_rhs, method="Radau")

    def solve_scipy_BDF():
        return _solve_scipy(rhs, method="BDF")

    def solve_scipy_numba_BDF():
        return _solve_scipy(numba_rhs, method="BDF")

    def solve_scipy_LSODA():
        return _solve_scipy(rhs, method="LSODA")

    def solve_scipy_numba_LSODA():
        return _solve_scipy(numba_rhs, method="LSODA")

    def solve_goss_ExplicitEuler():
        return goss_ExplicitEuler.solve(t)

    def solve_goss_RK2():
        return goss_RK2.solve(t)

    def solve_goss_RK4():
        return goss_RK4.solve(t)

    def solve_goss_RL1():
        return goss_RL1.solve(t)

    def solve_goss_RL2():
        return goss_RL2.solve(t)

    def solve_goss_GRL1():
        return goss_GRL1.solve(t)

    def solve_goss_GRL2():
        return goss_GRL2.solve(t)

    def solve_goss_RKF32():
        return goss_RKF32.solve(t)

    def solve_goss_ImplicitEuler():
        return goss_ImplicitEuler.solve(t)

    def solve_goss_ThetaSolver():
        return goss_ThetaSolver.solve(t)

    def solve_goss_ESDIRK23a():
        return goss_ESDIRK23a.solve(t)

    all_methods = {
        "goss (ExplicitEuler)": solve_goss_ExplicitEuler,
        "goss (RK2)": solve_goss_RK2,
        "goss (RK4)": solve_goss_RK4,
        "goss (RL1)": solve_goss_RL1,
        "goss (RL2)": solve_goss_RL2,
        "goss (GRL1)": solve_goss_GRL1,
        "goss (GRL2)": solve_goss_GRL2,
        "goss (ImplicitEuler)": solve_goss_ImplicitEuler,
        "goss (ThetaSolver)": solve_goss_ThetaSolver,
        "goss (RKF32)": solve_goss_RKF32,
        "goss (ESDIRK23a)": solve_goss_ESDIRK23a,
        "scipy (RK45)": solve_scipy_RK45,
        "scipy+numba (RK45)": solve_scipy_numba_RK45,
        "scipy (RK23)": solve_scipy_RK23,
        "scipy+numba (RK23)": solve_scipy_numba_RK23,
        "scipy (DOP853)": solve_scipy_DOP853,
        "scipy+numba (DOP853)": solve_scipy_numba_DOP853,
        "scipy (Radau)": solve_scipy_Radau,
        "scipy+numba (Radau)": solve_scipy_numba_Radau,
        "scipy (BDF)": solve_scipy_BDF,
        "scipy+numba (BDF)": solve_scipy_numba_BDF,
        "scipy (LSODA)": solve_scipy_LSODA,
        "scipy+numba (LSODA)": solve_scipy_numba_LSODA,
    }

    if selected_methods is None:
        selected_methods = list(all_methods.keys())

    def time_methods():

        time_outfile = Path(f"timings_{name}_{number}_{repeat}.npy")

        if recompute or not time_outfile.is_file():
            print("Time methods")
            timings = {}

            for label, solver in all_methods.items():
                if label not in selected_methods:
                    continue
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

        if recompute or not solutions_outfile.is_file():
            data = {}
            for label, solver in all_methods.items():
                if label not in selected_methods:
                    continue
                print(f"Solve {label}")
                y = solver()
                if y.shape[1] == len(t):
                    y = y.T

                data[label] = y
            np.save(solutions_outfile, data)

        data = np.load(solutions_outfile, allow_pickle=True).item()

        fig, ax = plt.subplots(len(selected_states), 1, figsize=(10, 10), sharex=True)
        if len(selected_states) == 1:
            ax = np.array([ax])
        lines = []
        labels = []
        for i, (label, y) in enumerate(data.items()):
            for idx, state_idx in enumerate(selected_states):
                (l,) = ax[idx].plot(t, y[:, state_idx])
                if i == 0:
                    ax[idx].set_ylabel(ode.state_names[state_idx])
            lines.append(l)
            labels.append(label)

        ax[len(selected_states) - 1].set_xlabel("$t$")

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
