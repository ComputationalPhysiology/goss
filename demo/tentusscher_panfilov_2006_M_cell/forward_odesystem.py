# type: ignore
import time

import goss
import matplotlib.pyplot as plt
import numpy as np
from gotran import load_ode


field_states = ["V", "Ca_i"]
field_parameters = ["g_CaL", "g_Kr"]
ode = goss.ParameterizedODE(
    load_ode("tentusscher_panfilov_2006_M_cell.oder"),
    field_states=field_states,
    field_parameters=field_parameters,
)

solver = goss.solvers.ESDIRK23a()

# Let us solve 9 ODEs at the same time
num_nodes = 9
system = goss.ODESystemSolver(num_nodes, solver, ode)
# Use 9 threads
system.num_threads = num_nodes

# Let us pick three block factors (1 representing baseline)
block_factors = np.array([1, 0.8, 0.6])
# And update the field parameters
field_parameters = system.field_parameters
# Block only CaL
field_parameters[:3, 0] *= block_factors
# Block only Kr
field_parameters[3:6, 1] *= block_factors
# Block both
field_parameters[6:9, 0] *= block_factors
field_parameters[6:9, 1] *= block_factors

# Update the field parameters
system.field_parameters = field_parameters

tstop = 5000.0
dt = 1.0
time_stamps = np.arange(0, tstop, dt)

field_states = np.zeros((time_stamps.size, num_nodes, 2))
field_states[0, :, :] = system.field_states

t0 = time.perf_counter()
for step, t in enumerate(time_stamps):
    system.forward(t, dt)
    field_states[step, :, :] = system.field_states
print(f"Elapsed time {time.perf_counter() - t0:.2f} seconds")

fig, ax = plt.subplots(2, 3, figsize=(8, 8), sharex=True, sharey="row")
for i, color in enumerate(["g", "b", "r"]):
    for j, linestyle in enumerate(["-", ":", "--"]):
        ax[0, i].plot(
            time_stamps,
            field_states[:, 3 * i + j, 0],
            linestyle=linestyle,
            label=f"block: {round((1-block_factors[j])*100)}%",
        )
        ax[1, i].plot(
            time_stamps,
            field_states[:, 3 * i + j, 1],
            linestyle=linestyle,
            label=f"block: {round((1-block_factors[j])*100)}%",
        )

for axi in ax.flatten():
    axi.legend()
    axi.grid()

ax[0, 0].set_title("CaL block")
ax[0, 1].set_title("Kr block")
ax[0, 2].set_title("CaL + Kr block")
ax[0, 0].set_ylabel("V")
ax[1, 0].set_ylabel("Cai")

plt.show()
