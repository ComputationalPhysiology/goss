import matplotlib.pyplot as plt
from gotran import load_ode
import goss
import numpy as np

oscilator = goss.ODE(load_ode("oscilator"))


solver = goss.solvers.ExplicitEuler(oscilator)
# solver = goss.solvers.RL1(oscilator)

tstop = 10.0
n_steps = 1000
time = np.linspace(0, tstop, n_steps + 1)

dt = time[1] - time[0]
u0 = np.array([1.0, 0.0])

u_exact = np.array([np.cos(time), np.sin(time)])
u1 = np.zeros_like(u_exact)
u2 = np.zeros_like(u_exact)
u1[:, 0] = u0
u2[:, 0] = u0

print(u1)
u = u0
for step in range(1, n_steps + 1):
    t = time[step]
    solver.forward(u, t, dt)
    u1[:, step] = u


fig, ax = plt.subplots()
ax.plot(time, u1[0, :], label="x")
ax.plot(time, u1[1, :], label="y")
plt.show()
