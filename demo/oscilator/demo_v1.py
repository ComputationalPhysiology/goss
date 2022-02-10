import goss
import matplotlib.pyplot as plt
import numpy as np
from gotran import load_ode

oscilator = goss.ODE(load_ode("oscilator"))


solver = goss.solvers.ExplicitEuler(oscilator)
# solver.set_parameter("ldt", 0.01)
# solver = goss.solvers.RL1(oscilator)

tstop = 10.0
n_steps = 1000
time = np.linspace(0, tstop, n_steps + 1)
dt = time[1] - time[0]


u_exact = np.array([np.cos(time), np.sin(time)])

u0 = np.array([1.0, 0.0])
u = np.zeros_like(u_exact)
u[:, 0] = u0

ui = u0
for step in range(1, n_steps + 1):
    t = time[step]
    solver.forward(ui, t, dt)
    u[:, step] = ui

fig, ax = plt.subplots()
ax.plot(time, u[0, :], label="u1 (computed)")
ax.plot(time, u[1, :], label="u2 (computed)")
ax.plot(time, u_exact[0, :], label="u1 (exact)", linestyle="--")
ax.plot(time, u_exact[1, :], label="u2 (exact)", linestyle="--")
ax.legend()
plt.show()
