# # Oscilator - Using generated code
#
# It might be beneficial to generate the C++ code first using the `gotran2goss` script in the CLI, and then use that code as input to the model. One advantage of this is that you can edit the C++ code directly before passing it to `goss`. This might be desirable if you want to implement logic into the model that is difficult to describe using `gotran` or if you want to bring in some third party library
#
# As before we first need to import the necessary libraries
#

import goss
import matplotlib.pyplot as plt
import numpy as np

#
# We assume that we ran the script
#
# ```
# goss gotran2goss oscilator.ode
# ```
#
# which produced the file `oscilator.h`. Inside this file there will be a method looking something like this
#
# ```c++
# extern "C" DLL_EXPORT goss::ParameterizedODE * create_Oscilator()
# {
#   return new goss::Oscilator;
# }
# ```
# This is the method that will create the class that is exported to the python interface. Note that the method is called `create_Oscilator`, but if you try to use a different ode-file the name might be different, but it will always be `create_{modelname}`. This name needs to be provided as a keyword argument to the constructor, so that in our case we will have
#
#

oscilator = goss.ParameterizedODE("oscilator.h", name="Oscilator")

#
# Next we create a solver
#

solver = goss.solvers.ThetaSolver(oscilator)

#
# solve the system
#

u0 = oscilator.get_ic()
dt = 0.01
t = np.arange(0, 10 + dt, dt)
u = solver.solve(t, y0=u0)
u_exact = np.array([np.cos(t), np.sin(t)])

#
# an finally plot the solution
#

fig, ax = plt.subplots()
ax.plot(t, u[:, 0], label="u1 (computed)")
ax.plot(t, u[:, 1], label="u2 (computed)")
ax.plot(t, u_exact[0, :], label="u1 (exact)", linestyle="--")
ax.plot(t, u_exact[1, :], label="u2 (exact)", linestyle="--")
ax.legend()
plt.show()

#
# ```{figure} _static/oscilator.png
# ---
# name: oscilator_v3
# ---
# Computed and exact solution of the oscilator ODE.
# ```
