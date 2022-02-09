from beatadjoint import *
from gotran import load_ode
from goss import dolfin_jit
from beatadjoint.gossplittingsolver import GOSSplittingSolver

# Parameters
Dx = Constant(1.5)
Dy = Constant(3.0)
D = as_tensor([[Dx, 0.], [0,Dy]])

dt = 0.125
tstop = 25.0
a = Constant(1.)
V_init = -85.
V_amp = 85.
t = Constant(0.)

# Domain and solution space
do_plot = False
L = 100.
N = 1024
#N = 128
domain = RectangleMesh(-L, -L, L, L, N, N)
cellmodel = dolfin_jit(load_ode("tentusscher_panfilov_2006_M_cell.ode"), field_states=["V"])
heart = CardiacModel(domain, t, D, None, cellmodel)
ps = GOSSplittingSolver.default_parameters()
ps["pde_solver"] = "monodomain"

ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
ps["MonodomainSolver"]["theta"] = 1.0
ps["ode_solver"]["solver"] = "RL1"
ps["ode_solver"]["num_threads"] = 8 / MPI.size(domain.mpi_comm())

# If cuda
ps["ode_solver"]["use_cuda"] = True
ps["ode_solver"]["cuda_params"]["float_precision"] = "double"
ps["ode_solver"]["cuda_params"]["solver"] = "rush_larsen"

ps["apply_stimulus_current_to_pde"] = True
solver = GOSSplittingSolver(heart, ps)

# Get solution fields
(u, um) = solver.solution_fields()

# Initialize function with different functions
# X planar
init_expr = Expression("V_amp/2.*(1-tanh(sqrt(a/(8*D))*(x[0]-x0)))+V_init", \
                       a=a, D=Dx, x0=-0.95*L, V_init=V_init, V_amp=V_amp)

init_expr = Expression("V_amp/2.*(1+tanh(-sqrt(a/(8*D))*2*sqrt(x[0]*x[0]+x[1]*x[1])/r))+V_init", r = 20., a=a, D=Dx, V_init=V_init, V_amp=V_amp)

u.interpolate(init_expr)

# Plot solution
if do_plot:
    plot(u, scale=1., range_min=-90., range_max=40.)

# Iterate
t = 0.
for timestep, (u, vm) in solver.solve((0, tstop), dt):

    # plot solution
    if do_plot:
        plot(u, scale=1., range_min=-90., range_max=40.)

    if MPI.rank(domain.mpi_comm()) == 0:
        print(timestep[0], "Min:", u.vector().min(), "Max:", u.vector().max())
    else:
        u.vector().min()
        u.vector().max()

if MPI.rank(domain.mpi_comm()) == 0:
    list_timings()

interactive()
