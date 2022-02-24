import goss
import matplotlib.pyplot as plt
import numpy as np
from gotran import load_ode

oscilator = goss.ODE(load_ode("oscilator"))


# solver = goss.solvers.ExplicitEuler(oscilator)
# solver = goss.solvers.RL1(oscilator)
solver = goss.solvers.ThetaSolver(oscilator)


u0 = oscilator.get_ic()
dt = 0.01
t = np.arange(0, 10 + dt, dt)
u = solver.solve(t, y0=u0)
u_exact = np.array([np.cos(t), np.sin(t)])

fig, ax = plt.subplots()
ax.plot(t, u[:, 0], label="u1 (computed)")
ax.plot(t, u[:, 1], label="u2 (computed)")
ax.plot(t, u_exact[0, :], label="u1 (exact)", linestyle="--")
ax.plot(t, u_exact[1, :], label="u2 (exact)", linestyle="--")
ax.legend()
plt.show()
