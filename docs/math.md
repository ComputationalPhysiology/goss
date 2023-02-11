# Mathematical background

If you the state of a process today and you know the physical laws that describes how this process changes, then you can formulate a differential equation to compute the state that that process a later time. For example if you know the whether today and the equations for modeling the change in the whether then we can predict the whether tomorrow.

Ordinary differential equations (ODEs) are equations of the form

```{math}
\frac{d \mathbf{x}}{d t} = f(\mathbf{x}, t; \mathbf{p})
```

Here $\mathbf{x}$ is a vector of state variables, $t$ is the time, $p$ are some known parameters, and $f$ describes how $\mathbf{x}$ evolves over time.

## A simple example
A very simple ODE is the model of an exponential growth

```{math}
\frac{{\rm d}u}{{\rm d}t} = a u
```
where $a$ is the growth constant. Here $u$ could for example be a population of some species that grows exponentially with a factor $a$, i.e at any point in time the change in $u$ is proportionally $u$ times some growth factor.

One way to solve this equations is to use the definition of the derivative
```{math}
\frac{{\rm d}u}{{\rm d}t} = \mathrm{lim}_{h \mapsto 0 }\frac{u(t + h) - u(t)}{h}
```
to find an approximate solution
```{math}
\frac{{\rm d}u}{{\rm d}t} \approx \frac{u(t + h) - u(t)}{h}
```
for a sufficiently small $h$. We can then instead write
```{math}
\frac{{\rm d}u}{{\rm d}t} = \approx \frac{u(t + h) - u(t)}{h} =  a u(t) \\
\implies u(t + h) = u(t) + h a u(t)
```
If we know the solution at time $t$ we can now approximate the solution at time $t + h$. This method is known as the _explicit Euler method_ of _forward Euler_. If we instead had evaluated the right hand side at $t+u$, i.e
```{math}
\frac{{\rm d}u}{{\rm d}t} = \approx \frac{u(t + h) - u(t)}{h} =  a u(t + h) \\
\implies u(t + h) = \frac{u(t)}{1 - ha}
```
we would obtain a different scheme known as the _implicit euler method_ or _backward Euler_.

Different schemes are useful in different situations. For examples for stiff ODEs, you might need to use an implicit scheme rather than an explicit scheme. In other cases it might be beneficial to apply adaptive time stepping to speed up the solving time.

## Schemes in `goss`
In `goss` we have implemented 11 different schemes which could


|          Name | Explicit/Implicit | Adaptive/Nonadaptive |
|---------------|-------------------|----------------------|
| ExplicitEuler | Explicit          | Nonadaptive          |
|           RK2 | Explicit          | Nonadaptive          |
|           RK4 | Explicit          | Nonadaptive          |
|           RL1 | Explicit          | Nonadaptive          |
|           RL2 | Explicit          | Nonadaptive          |
|          GRL1 | Explicit          | Nonadaptive          |
|          GRL2 | Explicit          | Nonadaptive          |
| ImplicitEuler | Implicit          | Nonadaptive          |
|   ThetaSolver | Implicit          | Nonadaptive          |
|         RKF32 | Explicit          | Adaptive             |
|     ESDIRK23a | Implicit          | Adaptive             |
