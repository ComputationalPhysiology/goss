# Comparison with Scipy and Scipy + Numba

In this section we have compared `goss` with [`scipy.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). We have also tested the cases when we use `numba` to jit compile the right hand side, which should speed up the solve time of the `scipy` solver if the right hand side is sufficiently complex (i.e has a lot of operations).

The setup for this benchmark is found in at [benchmarks/benchmark.py](https://github.com/ComputationalPhysiology/goss/blob/main/benchmarks/benchmark.py). In this benchmark we also perform some checks to make sure that the solution is properly converged so that we use the correct time steps for this different schemes.

The code for the right hand side using in the scipy solver was generated from the `.ode` file using `gotran`.

## Lorentz attractor

First benchmark is using the Lorentz attractor model, which is a model with three state variables

```{math}
\begin{align*}
\frac{dx}{dt} &= \sigma (y - x) \\
\frac{dy}{dt} &= x (\rho  - z) - y \\
\frac{dz}{dt} &= x y - \beta z
\end{align*}
```

with
- $\mathbf{x} = (x, y, z)$
- $\mathbf{p} = (\sigma, \rho, \beta)$


```{figure} _static/timings_lorentz_5_10.png
---
name: timings_tentusscher_panfilov_2006_M_cell_5_5
---

Timings for running the Lorentz attractor model running for 500 s and averaged over 10 runs
```




#  Ten Tusscher model

The next test case is the Ten Tusscher model {cite}`ten2006alternans`  which has 19 states

```{figure} _static/timings_tentusscher_panfilov_2006_M_cell_5_5.png
---
name: timings_tentusscher_panfilov_2006_M_cell_5_5
---

Timings for running the Ten Tusscher model running for 5000 ms and averaged over 5 runs
```


## References

```{bibliography}
:filter: docname in docnames
```
