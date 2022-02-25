import goss
import matplotlib.pyplot as plt
import numpy as np
from gotran import load_ode

# Import the ode in gotran
lorentz = load_ode("lorentz.ode")
# Jit compile the code for the right hand side
ode = goss.ParameterizedODE(lorentz)
# Select a solver and instantiate the solver
solver = goss.solvers.RKF32(ode)
# Select the time steps you want to solve for
t = np.linspace(0, 100, 10001)
# Solve the system
u = solver.solve(t)

# Plot the solution
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

ax.plot(u[:, 0], u[:, 1], u[:, 2], lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Lorenz Attractor")
fig.savefig("lorentz.png")
plt.show()
