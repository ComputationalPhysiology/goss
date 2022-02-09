from dolfin import *
from gotran import load_ode
from goss.dolfinutils import *

def setup_model(cellmodel_strs, domain, space):

    L = domain.coordinates().max()
    family, degree = family_and_degree_from_str(space)
    V = FunctionSpace(domain, family, degree)

    Vs = FunctionSpace(domain, "P", 1)

    field_parameters = dict(rice_model_2008=["hf", "hb"],
                            hybrid=["TMon_coop", "TMon_pow"])
    field_states = dict(rice_model_2008=["active"],
                        hybrid=["active"])#, "TCa"

    # Alter spatially varying paramters:
    param_scale = Expression("offset+scale*exp(-((x[0]-center_x)*(x[0]-center_x)+"\
                               "(x[1]-center_y)*(x[1]-center_y))/(sigma*sigma))",
                               center_x=3*L/4, center_y=L/4, offset=0.0, sigma=L/2, \
                               scale=1.0)

    cellmodels = [dolfin_jit(\
        load_ode(cellmodel), field_states=field_states[cellmodel], \
        field_parameters=field_parameters[cellmodel]) for cellmodel in cellmodel_strs]

    if "rice_model_2008" in cellmodel_strs:
        max_value = 0.14
        index = cellmodel_strs.index("rice_model_2008")
        model = cellmodels[index]

        values = dict(hf=0.03, hb=0.06)
        y_shift = dict(hf=2.5, hb=3)
        for param in ["hf", "hb"]:
            p0 = model.get_parameter(param)
            param_scale.offset = p0
            param_scale.scale = -(p0-values[param])
            param_scale.center_y = y_shift[param]*L/4
            p_func = Function(V, name=param)
            p_func.interpolate(param_scale)
            model.set_parameter(param, p_func)

    if "hybrid" in cellmodel_strs:
        max_value=0.14
        index = cellmodel_strs.index("hybrid")
        model = cellmodels[index]
        values = dict(TMon_coop=0.5, TMon_pow=0.5)
        y_shift = dict(TMon_coop=1.5, TMon_pow=2.0)
        for param in list(values.keys()):
            p0 = model.get_parameter(param)
            param_scale.offset = p0
            param_scale.scale = -(p0-values[param])
            param_scale.center_y = y_shift[param]*L/4
            p_func = Function(V, name=param)
            p_func.interpolate(param_scale)
            model.set_parameter(param, p_func)

    if len(cellmodels)==2:
        subdomain = CompiledSubDomain("x[1] <= 0.5")
        if space == "P_1":
            cellmodel_domains = VertexFunction("size_t", domain, 10)
        else:
            cellmodel_domains = CellFunction("size_t", domain, 10)
        subdomain.mark(cellmodel_domains, 20)
        cellmodels = dict(label_model for label_model in zip([10,20], cellmodels))

    elif len(cellmodels)==1:
        cellmodels = cellmodels[0]
        cellmodel_domains=None

    else:
        assert False

    params = DOLFINODESystemSolver.default_parameters()
    params["solver"] = "RL1"
    solver = DOLFINODESystemSolver(domain, cellmodels, domains=cellmodel_domains, \
                                   space=space, params=params)
    u = Function(solver.state_space)
    solver.from_field_states(u)

    if space != "P1":
        u_plot = Function(Vs)
        if solver.num_field_states > 1:
            u_plot.assign(project(u[0], Vs))
        else:
            u_plot.assign(project(u, Vs))

    elif solver.num_field_states > 1:
        u_plot = u.split(True)[0]
    else:
        u_plot = u

    t=0
    dt=.1
    tstop=300

    while t < tstop:
        solver.step((t,t+dt), u)
        if (t%10.)<dt:
            print("t:", t)
            if space != "P1":
                if solver.num_field_states > 1:
                    u_plot.assign(project(u[0], Vs))
                else:
                    u_plot.assign(project(u, Vs))

            elif solver.num_field_states > 1:
                assign(u_plot, u.sub(0))

            plot(u_plot, scale=0., range_max=max_value, range_min=0.)

        t+=dt

    plot(u_plot, scale=0., range_max=max_value, range_min=0., interactive=True)


if __name__ == "__main__":

    cellmodel_strs = ["rice_model_2008", "hybrid"]
    space = "Quadrature_2"#"DG_0"#, "Quadrature_2"
    space = "P_1"#"DG_0"#, "Quadrature_2"

    # Define mesh
    domain = UnitSquareMesh(10, 10)
    setup_model(cellmodel_strs, domain, space)
