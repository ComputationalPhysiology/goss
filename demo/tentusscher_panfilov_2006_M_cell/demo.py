import goss
import matplotlib.pyplot as plt
from gotran import load_ode

gotran_ode = load_ode("tentusscher_panfilov_2006_M_cell.ode")
ode = goss.ODE(gotran_ode)


# solver = goss.solvers.ExplicitEuler(ode)
# solver = goss.solvers.RL1(ode)
# solver = goss.solvers.GRL1(ode)
solver = goss.solvers.ESDIRK23a(ode)


y, t = solver.solve(0, 1000, dt=1.0)

V_index = gotran_ode.state_symbols.index("V")
Cai_index = gotran_ode.state_symbols.index("Ca_i")

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(t, y[:, V_index])
ax[1].plot(t, y[:, Cai_index])
ax[0].set_title("V")
ax[1].set_title("Cai")
plt.show()
