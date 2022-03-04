import dolfin
import cbcbeat
import math
import logging
import goss
import tqdm
import gotran
from cbcbeat.gossplittingsolver import GOSSplittingSolver

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


domain = dolfin.UnitSquareMesh(100, 100)
domain.coordinates()[:] *= 10


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


gotran_ode = gotran.load_ode("tentusscher_panfilov_2006_M_cell.ode")


stim_amplitude = 50.0
stim_duration = 1.0

dt = 0.5
T = 100

cellmodel = goss.dolfinutils.DOLFINParameterizedODE(
    gotran_ode,
    field_states=["V"],
)
# Do not apply any stimulus in the cell model
cellmodel.set_parameter("stim_amplitude", 0)

V = dolfin.FunctionSpace(domain, "CG", 1)
L = domain.coordinates().max()


chi = 12000.0  # cm^{-1}
s_il = 300.0 / chi  # mS
s_it = s_il / 2  # mS
s_el = 200.0 / chi  # mS
s_et = s_el / 1.2  # mS

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

heart = cbcbeat.CardiacModel(domain, time, M_i, M_e, cellmodel, stimulus)


solver = GOSSplittingSolver(heart, ps)

(vs_, vs, vur) = solver.solution_fields()
V = vs.function_space()
v = dolfin.Function(V)

# Solve
solutions = solver.solve((0, T), dt)

vfile = dolfin.XDMFFile(dolfin.MPI.comm_world, "v.xdmf")
for ((t0, t1), fields) in tqdm.tqdm(solutions, total=int(T / dt)):
    v.assign(vs_)
    vfile.write(v, t0)
