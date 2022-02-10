import goss
import matplotlib.pyplot as plt
import numpy as np
from gotran import load_ode

oscilator = goss.ODE(load_ode("oscilator"))


# solver = goss.solvers.ExplicitEuler(oscilator)
# solver = goss.solvers.RL1(oscilator)
solver = goss.solvers.ThetaSolver(oscilator)


u0 = oscilator.get_ic()
u, time = solver.solve(0, 10.0, dt=0.01, y0=u0)
u_exact = np.array([np.cos(time), np.sin(time)])

fig, ax = plt.subplots()
ax.plot(time, u[:, 0], label="u1 (computed)")
ax.plot(time, u[:, 1], label="u2 (computed)")
ax.plot(time, u_exact[0, :], label="u1 (exact)", linestyle="--")
ax.plot(time, u_exact[1, :], label="u2 (exact)", linestyle="--")
ax.legend()
plt.show()
