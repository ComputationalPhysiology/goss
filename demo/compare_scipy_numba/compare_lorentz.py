import benchmark
import numpy as np


def rhs(t, u, sigma, rho, beta):
    x, y, z = u
    return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]


def main():
    benchmark.main(
        rhs,
        "lorentz.ode",
        t=np.linspace(0, 100, num=1001),
        args=(10.0, 28.0, 8 / 3),
        run_plot=True,
        run_timings=False,
    )


if __name__ == "__main__":
    main()
