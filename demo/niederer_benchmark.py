# (bidomain)=
# # Bidomain model
#

import dolfin
import cbcbeat
import goss
import numpy as np
import gotran
from cbcbeat.gossplittingsolver import GOSSplittingSolver

ic = {
    "V": -85.23,  # mV
    "Xr1": 0.00621,
    "Xr2": 0.4712,
    "Xs": 0.0095,
    "m": 0.00172,
    "h": 0.7444,
    "j": 0.7045,
    "d": 3.373e-05,
    "f": 0.7888,
    "f2": 0.9755,
    "fCass": 0.9953,
    "s": 0.999998,
    "r": 2.42e-08,
    "Ca_i": 0.000126,  # millimolar
    "R_prime": 0.9073,
    "Ca_SR": 3.64,  # millimolar
    "Ca_ss": 0.00036,  # millimolar
    "Na_i": 8.604,  # millimolar
    "K_i": 136.89,  # millimolar
}


def setup_model(dx=0.5):

    Lx = 20.0  # mm
    Ly = 7.0  # mm
    Lz = 3.0  # mm

    N = lambda v: int(np.rint(v))
    mesh = dolfin.BoxMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0, 0.0),
        dolfin.Point(Lx, Ly, Lz),
        N(Lx / dx),
        N(Ly / dx),
        N(Lz / dx),
    )

    L = mesh.coordinates().max()

    # Surface to volume ratio
    chi = 140.0  # mm^{-1}
    # Membrane capacitance
    C_m = 0.01  # mu F / mm^2

    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = 0.17  # mS / mm
    sigma_it = 0.019  # mS / mm
    sigma_el = 0.62  # mS / mm
    sigma_et = 0.24  # mS / mm

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(C_m * chi)
    s_l = sigma_l / (C_m * chi)  # mm^2 / ms
    s_t = sigma_t / (C_m * chi)  # mm^2 / ms

    M = dolfin.as_tensor(((s_l, 0, 0), (0, s_t, 0), (0, 0, s_t)))

    # Define time
    time = dolfin.Constant(0.0)

    # Mark stimulation region defined as [0, L]^3
    S1_marker = 1
    L = 1.5
    S1_subdomain = dolfin.CompiledSubDomain(
        "x[0] <= L + DOLFIN_EPS && x[1] <= L + DOLFIN_EPS && x[2] <= L + DOLFIN_EPS",
        L=L,
    )
    S1_markers = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    S1_subdomain.mark(S1_markers, S1_marker)

    # Define stimulation (NB: region of interest carried by the mesh
    # and assumptions in cbcbeat)
    duration = 2.0  # ms
    A = 50000.0  # mu A/cm^3
    cm2mm = 10.0
    factor = 1.0 / (chi * C_m)  # NB: cbcbeat convention
    amplitude = factor * A * (1.0 / cm2mm) ** 3  # mV/ms
    I_s = dolfin.Expression(
        "time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=0.0,
        duration=duration,
        amplitude=amplitude,
        degree=0,
    )
    # Store input parameters in cardiac model
    stimulus = cbcbeat.Markerwise((I_s,), (1,), S1_markers)
    return mesh, time, M, stimulus


def run_benchmark_goss(mesh, time, M, stimulus, T, dt, scheme="GRL1"):
    gotran_ode = gotran.load_ode("tentusscher_panfilov_2006_M_cell.ode")

    cellmodel = goss.dolfinutils.DOLFINParameterizedODE(
        gotran_ode,
        field_states=gotran_ode.state_symbols,
    )
    # Do not apply any stimulus in the cell model
    cellmodel.set_parameter("stim_amplitude", 0)
    cellmodel.set_initial_conditions(**ic)

    heart = cbcbeat.CardiacModel(mesh, time, M, None, cellmodel, stimulus)

    ps = GOSSplittingSolver.default_parameters()
    ps["pde_solver"] = "monodomain"
    ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
    ps["MonodomainSolver"]["preconditioner"] = "sor"
    ps["MonodomainSolver"]["theta"] = 0.5
    ps["MonodomainSolver"]["default_timestep"] = dt
    ps["MonodomainSolver"]["use_custom_preconditioner"] = False

    ps["theta"] = 0.5
    ps["enable_adjoint"] = False
    ps["apply_stimulus_current_to_pde"] = True
    ps["ode_solver"]["solver"] = scheme

    V_index = heart.cell_models().state_names.index("V")

    solver = GOSSplittingSolver(heart, ps, V_index=V_index)

    # Extract the solution fields and set the initial conditions
    (vs_, vs, vur) = solver.solution_fields()

    # Set-up separate potential function for post processing
    VS0 = vs.function_space().sub(V_index)
    V = VS0.collapse()
    v = dolfin.Function(V)

    # Set-up object to optimize assignment from a function to subfunction
    assigner = dolfin.FunctionAssigner(V, VS0)
    assigner.assign(v, vs_.sub(V_index))

    # Solve
    solutions = solver.solve((0, T), dt)

    vfile = dolfin.XDMFFile(dolfin.MPI.comm_world, "niederer_benchmark_goss.xdmf")
    vfile.write(v, 0.0)
    for (i, ((t0, t1), fields)) in enumerate(solutions):
        if (i % 20 == 0) and dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
            print("Reached t=%g/%g, dt=%g" % (t0, T, dt))
        assigner.assign(v, vs_.sub(V_index))
        vfile.write(v, t1)


def run_benchmark_vanilla(mesh, time, M, stimulus, T, dt, scheme="GRL1"):

    CellModel = cbcbeat.Tentusscher_panfilov_2006_epi_cell

    ps = cbcbeat.SplittingSolver.default_parameters()
    ps["pde_solver"] = "monodomain"
    ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
    ps["MonodomainSolver"]["theta"] = 0.5
    ps["MonodomainSolver"]["preconditioner"] = "sor"
    ps["MonodomainSolver"]["default_timestep"] = dt
    ps["MonodomainSolver"]["use_custom_preconditioner"] = False
    ps["theta"] = 0.5
    ps["apply_stimulus_current_to_pde"] = True
    ps["CardiacODESolver"]["scheme"] = scheme

    cellmodel = CellModel(init_conditions=ic)

    # Set-up cardiac model
    heart = cbcbeat.CardiacModel(mesh, time, M, None, cellmodel, stimulus)

    solver = cbcbeat.SplittingSolver(heart, ps)

    # Extract the solution fields and set the initial conditions
    (vs_, vs, vur) = solver.solution_fields()
    vs_.assign(cellmodel.initial_conditions())

    # Set-up separate potential function for post processing
    VS0 = vs.function_space().sub(0)
    V = VS0.collapse()
    v = dolfin.Function(V)

    # Set-up object to optimize assignment from a function to subfunction
    assigner = dolfin.FunctionAssigner(V, VS0)
    assigner.assign(v, vs_.sub(0))

    t0 = 0.0
    solutions = solver.solve((t0, T), dt)
    vfile = dolfin.XDMFFile(dolfin.MPI.comm_world, "niederer_benchmark_vanilla.xdmf")
    vfile.write(v, 0.0)
    for (i, ((t0, t1), fields)) in enumerate(solutions):
        if (i % 20 == 0) and dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
            print("Reached t=%g/%g, dt=%g" % (t0, T, dt))
        assigner.assign(v, vs_.sub(0))
        vfile.write(v, t1)


def main():

    dt = 0.5
    dx = 0.5
    T = 100.0
    mesh, time, M, stimulus = setup_model(dx=dx)

    with dolfin.Timer("Goss"):
        run_benchmark_goss(mesh=mesh, time=time, M=M, stimulus=stimulus, dt=dt, T=T)

    with dolfin.Timer("Vanilla"):
        run_benchmark_vanilla(mesh=mesh, time=time, M=M, stimulus=stimulus, dt=dt, T=T)

    dolfin.list_timings(dolfin.TimingClear.keep, [dolfin.TimingType.wall])


if __name__ == "__main__":
    main()
