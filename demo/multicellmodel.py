"""
This test tests the splitting solver for the bidomain equations with a
FitzHughNagumo and ten-Tusscher model.
"""

import math

import dolfin
import cbcbeat

from goss.dolfinutils import DOLFINParameterizedODE

from cbcbeat.gossplittingsolver import GOSSplittingSolver

dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["optimize"] = True
dolfin.parameters["form_compiler"]["quadrature_degree"] = 2


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


# Set-up solver
ps = GOSSplittingSolver.default_parameters()
ps["pde_solver"] = "monodomain"

ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
ps["MonodomainSolver"]["theta"] = 0.5

ps["theta"] = 0.5
ps["enable_adjoint"] = False
ps["apply_stimulus_current_to_pde"] = True
ps["ode_solver"]["solver"] = "RL1"
ps["ode_solver"]["num_threads"] = 0

harmonic_mean = bool(ps["pde_solver"] == "monodomain")

domain = dolfin.UnitSquareMesh(100, 100)
domain.coordinates()[:] *= 10
membrane_potential = "V"


stim_amplitude = 50.0
stim_duration = 1.0

dt = 0.5
T = 100.0

labels = dict(tentusscher_panfilov_2006_epi_cell=10, fitzhughnagumo=20)

field_parameters = dict(
    tentusscher_panfilov_2006_epi_cell=["g_CaL"],
    fitzhughnagumo=["b"],
)
field_states = dict(tentusscher_panfilov_2006_epi_cell=["V"], fitzhughnagumo=["V"])

cellmodels = [
    DOLFINParameterizedODE(
        cellmodel,
        field_states=field_states[cellmodel],
        field_parameters=field_parameters[cellmodel],
    )
    for cellmodel in labels.keys()
]

for cellmodel in cellmodels:
    cellmodel.set_initial_conditions(V=-86.2)

# Create subdomains for applying the different ODEs
subdomain = dolfin.CompiledSubDomain("x[0] <= 5.0")
cellmodel_domains = dolfin.MeshFunction("size_t", domain, 0, 10)
subdomain.mark(cellmodel_domains, 20)
cellmodels = cbcbeat.MultiCellModel(
    cellmodels,
    list(labels.values()),
    cellmodel_domains,
)

# Create scalar FunctionSpace
V = dolfin.FunctionSpace(domain, "CG", 1)
L = domain.coordinates().max()

# Alter spatially varying paramters:
param_scale = dolfin.Expression(
    "offset+scale*exp(-((x[0]-center_x)*(x[0]-center_x)+"
    "(x[1]-center_y)*(x[1]-center_y))/(sigma*sigma))",
    center_x=3 * L / 4,
    center_y=L / 4,
    offset=0.0,
    sigma=L / 2,
    scale=1.0,
    degree=2,
)

# Tentusscher parameter
ode_tt = cellmodels[labels["tentusscher_panfilov_2006_epi_cell"]]
g_CaL_0 = ode_tt.get_parameter("g_CaL")
param_scale.offset = g_CaL_0
param_scale.scale = -g_CaL_0 * 0.95
param_scale.center_y = 3 * L / 4
g_CaL_func = dolfin.Function(V)
g_CaL_func.interpolate(param_scale)
ode_tt.set_parameter("g_CaL", g_CaL_func)

# FHN paramter
ode_fhn = cellmodels[labels["fitzhughnagumo"]]

# Set-up cardiac model
k = 0.00004
V_rest = -85.0
V_threshold = -70.0
V_peak = 40.0
V_amp = V_peak - V_rest
l = 0.63  # noqa: E741
b = 0.013
param_scale.scale = b
param_scale.offset = b
param_scale.center_x = L / 4
param_scale.center_y = L / 4
b_func = dolfin.Function(V)
b_func.interpolate(param_scale)
cell_parameters = {
    "c_1": k * V_amp**2,
    "c_2": k * V_amp,
    "c_3": b / l,
    "a": (V_threshold - V_rest) / V_amp,
    "b": b_func,
    "V_rest": V_rest,
    "V_peak": V_peak,
}
# Set FHN specific parameters
for params in cell_parameters.items():
    ode_fhn.set_parameter(*params)

# Define conductivities
chi = 12000.0  # cm^{-1}
s_il = 300.0 / chi  # mS
s_it = s_il / 2  # mS
s_el = 200.0 / chi  # mS
s_et = s_el / 1.2  # mS

# Conductivities
if harmonic_mean:
    sl = s_il * s_el / (s_il + s_el)
    st = s_it * s_et / (s_it + s_et)
    M_i = dolfin.as_tensor(((dolfin.Constant(sl), 0), (0, dolfin.Constant(st))))
else:
    M_i = dolfin.as_tensor(((dolfin.Constant(s_il), 0), (0, dolfin.Constant(s_it))))

M_e = dolfin.as_tensor(((dolfin.Constant(s_el), 0), (0, dolfin.Constant(s_et))))

# Stimulus
stim_marker = 1
domain_size = domain.coordinates().max()
stim_subdomain = StimSubDomain((domain_size / 2.0, 0.0), domain_size / 5.0)
stim_domain = dolfin.MeshFunction("size_t", domain, domain.topology().dim(), 0)
stim_subdomain.mark(stim_domain, stim_marker)
time = dolfin.Constant(0.0)
stim = dolfin.Expression(
    "time > start ? (time <= (duration + start) ? " "amplitude : 0.0) : 0.0",
    time=time,
    duration=stim_duration,
    start=1.0,
    amplitude=stim_amplitude,
    degree=2,
)
stimulus = cbcbeat.Markerwise([stim], [stim_marker], stim_domain)

# Create and return CardiacModel
heart = cbcbeat.CardiacModel(domain, time, M_i, M_e, cellmodels, stimulus)

solver = GOSSplittingSolver(heart, ps)

# Extract the solution fields and set the initial conditions
(vs_, vs, vur) = solver.solution_fields()

# Set-up separate potential function for post processing
V = vs.function_space()
v = dolfin.Function(V)

total = dolfin.Timer("Total solver time")
solutions = solver.solve((0, T), dt)

vfile = dolfin.XDMFFile(dolfin.MPI.comm_world, "multicellmodel.xdmf")
vfile.write(v, 0.0)
for (i, ((t0, t1), fields)) in enumerate(solutions):
    if (i % 20 == 0) and dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
        print("Reached t=%g/%g, dt=%g" % (t0, T, dt))
    v.assign(vs_)
    vfile.write(v, t1)
