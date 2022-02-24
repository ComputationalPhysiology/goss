import benchmark
import numpy as np


def rhs(t, u, sigma, rho, beta):
    x, y, z = u
    return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]


selected_methods = [
    "goss (ExplicitEuler)",
    "goss (RK2)",
    "goss (RK4)",
    "goss (RL1)",
    "goss (RL2)",
    "goss (GRL1)",
    "goss (GRL2)",
    "goss (ImplicitEuler)",
    "goss (ThetaSolver)",
    "goss (RKF32)",
    # "goss (ESDIRK23a)",
    "scipy (RK45)",
    "scipy+numba (RK45)",
    "scipy (RK23)",
    "scipy+numba (RK23)",
    "scipy (DOP853)",
    "scipy+numba (DOP853)",
    "scipy (Radau)",
    "scipy+numba (Radau)",
    "scipy (BDF)",
    "scipy+numba (BDF)",
    "scipy (LSODA)",
    "scipy+numba (LSODA)",
]


def main():
    benchmark.main(
        rhs,
        "lorentz.ode",
        t=np.linspace(0, 100, num=1001),
        args=(10.0, 28.0, 8 / 3),
        run_plot=True,
        run_timings=False,
        recompute=True,
        selected_methods=selected_methods,
        selected_states=[
            0,
            1,
        ],
        internal_time_step=0.001,
        number=1,
        repeat=2,
    )


if __name__ == "__main__":
    main()
