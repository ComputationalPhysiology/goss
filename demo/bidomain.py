# (bidomain)=
# # Bidomain model
#
# In this demo we solve the Bidomain  model
# ```{math}
#  \nabla \cdot \left( M_i \nabla v \right) + \nabla \cdot \left( M_i \nabla u_e \right)= \chi C_m \frac{\partial v}{\partial t} + \chi I_{ion} \\
#  \nabla \cdot \left( M_i \nabla v \right) + \nabla \cdot \left( (M_i + M_e) \nabla u_e \right) = 0
# ```
# on a square domain, $\Omega = [0, 10] \times [0, 10]$. Here subscript `i` means *intracellular* and subscript `e` means extracellular. This is true for conductivities ($M_i$ and $M_e$) and the potential ($u_i$ and $u_e$). We also have the membrane potential $v = u_i - u_e$ which is the potential difference in the intracellular and extracellular domains.
#
# Further we have the area of the cell membrane per unit volume $\chi$, the capacitance of the cell membrane $C_m$ and the Ionic current across the membrane.
#
# For a full derivation of the Bidomain equations see {cite}`sundnes2007computing` or {cite}`keener2009mathematical`.
#
# In this demo we will use the Tentusscher model {cite}`ten2006alternans` to model the ionic current across the membrane.
#
# First we need to make the necessary imports
#

import dolfin
import cbcbeat
import math
import logging
import goss
import tqdm
import gotran
from cbcbeat.gossplittingsolver import GOSSplittingSolver

# Note that we also import the `GOSSplittingSolver` from [`cbcbeat`](https://github.com/ComputationalPhysiology/cbcbeat). `cbcbeat` is a library for solving the problem in cardiac electrophysiology.

# First we initialize some parameters for the solver

ps = GOSSplittingSolver.default_parameters()
ps["pde_solver"] = "bidomain"
# ps["BidomainSolver"]["use_avg_u_constraint"] = False
# ps["BidomainSolver"]["linear_solver_type"] = "direct"
ps["BidomainSolver"]["linear_solver_type"] = "iterative"
ps["BidomainSolver"]["theta"] = 0.5
ps["theta"] = 0.5
ps["enable_adjoint"] = False
ps["apply_stimulus_current_to_pde"] = True
ps["ode_solver"]["solver"] = "GRL1"
ps["ode_solver"]["num_threads"] = 0
dolfin.set_log_level(logging.WARNING)

# We define the domain, i.e a unit square that we scale by a factor of 10.0

domain = dolfin.UnitSquareMesh(100, 100)
domain.coordinates()[:] *= 10


# Next we define the stimulus domain which will be a circle with


class StimSubDomain(dolfin.SubDomain):
    def __init__(self, center, radius):
        self.x0, self.y0 = center
        self.radius = radius
        super().__init__()

    def inside(self, x, on_boundary):
        r = math.sqrt((x[0] - self.x0) ** 2 + (x[1] - self.y0) ** 2)
        if r < self.radius:
            return True
        return False


# Create a mesh function containing markers for the stimulus domain
stim_domain = dolfin.MeshFunction("size_t", domain, domain.topology().dim(), 0)
# Set all markers to zero
stim_domain.set_all(0)
# We mark the stimulus domain with a different marker
stim_marker = 1
domain_size = domain.coordinates().max()
# Create a domain
stim_subdomain = StimSubDomain(
    center=(domain_size / 2.0, 0.0),
    radius=domain_size / 5.0,
)
# And mark the domain
stim_subdomain.mark(stim_domain, stim_marker)

# Next we create the stimulus protocol

# Strength of the amplitude (this is based on the ionic model)
stim_amplitude = 50.0
# Duration of the stimulus
stim_duration = 1.0
# Make a constant represting time
time = dolfin.Constant(0.0)
# Make an expression that we apply a stimulus at 1 ms for a duration of 1 ms
stim = dolfin.Expression(
    "time > start ? (time <= (duration + start) ? " "amplitude : 0.0) : 0.0",
    time=time,
    duration=stim_duration,
    start=1.0,
    amplitude=stim_amplitude,
    degree=2,
)
# Make a stimulus object to be passed to cbcbeat
stimulus = cbcbeat.Markerwise([stim], [stim_marker], stim_domain)


# We load the ode representing the ionic model

gotran_ode = gotran.load_ode("tentusscher_panfilov_2006_M_cell.ode")

# and create a cell model with the membrane potential as a field state
# This will make the membrane potential being spatial dependent-
cellmodel = goss.dolfinutils.DOLFINParameterizedODE(
    gotran_ode,
    field_states=["V"],
)

# Do not apply any stimulus in the cell model since we will do this through the bidomain model
cellmodel.set_parameter("stim_amplitude", 0)

# Define the conductivities (where does these values come from?)

chi = 12000.0  # cm^{-1}
s_il = 300.0 / chi  # mS
s_it = s_il / 2  # mS
s_el = 200.0 / chi  # mS
s_et = s_el / 1.2  # mS

# and the conductivity tensors
#
# ```{math}
# M_i = \begin{pmatrix}
# s_{il} & 0 \\
# 0 & s_{it}
# \end{pmatrix}, \;\;
# M_e = \begin{pmatrix}
# s_{el} & 0 \\
# 0 & s_{et}
# \end{pmatrix}
# ```

M_i = dolfin.as_tensor(((dolfin.Constant(s_il), 0), (0, dolfin.Constant(s_it))))
M_e = dolfin.as_tensor(((dolfin.Constant(s_el), 0), (0, dolfin.Constant(s_et))))


# Create and the CardiacModel in `cbcbeat`

heart = cbcbeat.CardiacModel(domain, time, M_i, M_e, cellmodel, stimulus)

# and initialize the solver
solver = GOSSplittingSolver(heart, ps)

# We extract the membrane potential from the solution fields

(vs_, vs, vur) = solver.solution_fields()
V = vs.function_space()
v = dolfin.Function(V)

# and solve the model for 100 ms with increments of 0.5. We also save the resulting membrane potential in and XDMF file that can be visualized in Paraview.

dt = 0.5
T = 100
vfile = dolfin.XDMFFile(dolfin.MPI.comm_world, "bidomain.xdmf")
for ((t0, t1), fields) in tqdm.tqdm(solver.solve((0, T), dt), total=int(T / dt)):
    v.assign(vs_)
    vfile.write(v, t0)


#
# ```{figure} _static/bidomain.png
# ---
# name: bidomain_fig
# ---
# Membrane potential
# ```
#


# ## Reference
#
# ```{bibliography}
# :filter: docname in docnames
# ```
