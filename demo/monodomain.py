# (monodomain)=
# # Monodomain model
#
# In this demo we solve the Monodomain model
# ```{math}
# \frac{\lambda}{1 + \lambda} \nabla \cdot \left( M_i \nabla v \right) = \chi C_m \frac{\partial v}{\partial t} + \chi I_{ion} \; x \in \Omega
# ```
# on a square domain, $\Omega = [0, 10] \times [0, 10]$. The monodomain model is a simplified version of the {ref}`bidomain` where we assume that the extracellular conductivity is proportional to the intracellular conductivity
# ```{math}
# M_e = \lambda M_i
# ```
# where $\lambda$ is a constant scalar.
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

# First we initialize some parameters for the solver

ps = GOSSplittingSolver.default_parameters()
ps["pde_solver"] = "monodomain"
ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
ps["MonodomainSolver"]["theta"] = 0.5
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
#


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

# Next we need to define the conductivity tensors. Here we use the same parameters as in {ref}`bidomain`

chi = 12000.0  # cm^{-1}
s_il = 300.0 / chi  # mS
s_it = s_il / 2  # mS
s_el = 200.0 / chi  # mS
s_et = s_el / 1.2  # mS

# and make a new conductivity tensor by taking the harmonic mean of the parameters

sl = s_il * s_el / (s_il + s_el)
st = s_it * s_et / (s_it + s_et)
M_i = dolfin.as_tensor(((dolfin.Constant(sl), 0), (0, dolfin.Constant(st))))


# Create and the CardiacModel in `cbcbeat`. Since we use a monodomain model we don't need to pass the extracellular conductivity

heart = cbcbeat.CardiacModel(
    domain=domain,
    time=time,
    M_i=M_i,
    M_e=None,
    cell_models=cellmodel,
    stimulus=stimulus,
)

# and initialize the solver
solver = GOSSplittingSolver(heart, ps)

# We extract the membrane potential from the solution fields

(vs_, vs, vur) = solver.solution_fields()
V = vs.function_space()
v = dolfin.Function(V)

# and solve the model for 100 ms with increments of 0.5. We also save the resulting membrane potential in and XDMF file that can be visualized in Paraview.

dt = 0.5
T = 100
vfile = dolfin.XDMFFile(dolfin.MPI.comm_world, "monodomain.xdmf")
for ((t0, t1), fields) in tqdm.tqdm(solver.solve((0, T), dt), total=int(T / dt)):
    v.assign(vs_)
    vfile.write(v, t0)

#
# ```{figure} _static/monodomain.png
# ---
# name: monodomain_fig
# ---
# Membrane potential
# ```
#


# ## Reference
#
# ```{bibliography}
# :filter: docname in docnames
# ```
#
