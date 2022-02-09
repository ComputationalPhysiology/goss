# Copyright (C) 2013 Johan Hake
#
# This file is part of GOSS.
#
# GOSS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GOSS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GOSS. If not, see <http://www.gnu.org/licenses/>.

__all__ = ["DOLFINODESystemSolver", "dolfin_jit", "family_and_degree_from_str"]

from collections import OrderedDict
import numpy as np
import types

# from . import cpp
from .compilemodule import jit as goss_jit

try:
    import dolfin as d
except ImportError:
    raise ImportError("dolfin is not present")

# Check version
from distutils.version import LooseVersion

if LooseVersion(d.__version__) <= LooseVersion("1.4.0"):
    raise ImportError("dolfin version need to be 1.4.0 or higher")

# Import Gotran and try import cuda solver
from goss import goss_solvers
from goss.solvers import ImplicitODESolver
import goss.cuda

enable_cuda = goss.cuda.cuda is not None


def entity_to_dofs(V):
    assert isinstance(V, d.FunctionSpaceBase)
    mesh = V.mesh()
    dim = mesh.topology().dim()
    dm = V.dofmap()
    if V.ufl_element().family() == "Lagrange":
        return d.vertex_to_dof_map(V)

    # FIXME: np.uintp?
    # Create an array
    num_entity_dofs = dm.num_entity_dofs(dim)
    entity_to_dof = np.zeros(mesh.num_cells() * num_entity_dofs, np.intc)
    num_fields = V.num_sub_spaces()
    if num_fields > 0:
        dms = [V.sub(i).dofmap() for i in range(num_fields)]
    else:
        num_fields = 1
        dms = [dm]

    for cell in d.cells(mesh):
        index = cell.index()
        for field_ind, dm_l in enumerate(dms):
            entity_to_dof[
                index * num_entity_dofs
                + field_ind : (index + 1) * num_entity_dofs
                + field_ind : num_fields
            ] = dms[field_ind].cell_dofs(index)

    return entity_to_dof


def _set_parameter(self, name, value):
    """
    Set the value of a parameter

    Argument:
    name : str
      The name of the paramter
    value : scalar, Function
      The value of the parameter. If the value is a Function the parameter
      needs to be a field parameter.
    """
    from . import _gosscpp

    assert isinstance(name, str), "expected a str as the name argument"
    assert isinstance(
        value,
        (int, float, d.Function),
    ), "expected a scalar or a Function for the value argument"

    # If the value is a scalar just call the original set_parameter
    if isinstance(value, (float, int)):
        _gosscpp.ParameterizedODE.set_parameter(self, name, value)
        self.changed_scalar_parameters[name] = value
    else:
        field_param_names = self.get_field_parameter_names()
        assert name in field_param_names, (
            "'%s' is not a field parameter in this ode" % name
        )
        self.field_params[name] = value
        self.changed_field_parameters.append(name)


def dolfin_jit(
    ode,
    field_states=None,
    field_parameters=None,
    monitored=None,
    code_params=None,
    cppargs=None,
):
    """
    Generate a goss::ODEParameterized from a gotran ode and JIT
    compile it. Add help methods to the jit compiled ODE to set field
    parameters from a dolfin Function

    Arguments:
    ----------
    ode : gotran.ODE
        The gotran ode, either as an ODE or as an ODERepresentation
    field_states : list
        A list of state names, which should be treated as field states
    field_parameters : list
        A list of parameter names, which should be treated as field parameters
    monitored : list
        A list of names of intermediates of the ODE. Code for monitoring
        the intermediates will be generated.
    code_params : dict
        Parameters controling the code generation
    cppargs : str
        Default C++ argument passed to the C++ compiler
    """
    from . import _gosscpp

    # Compile ode
    compiled_ode = goss_jit(
        ode,
        field_states,
        field_parameters,
        monitored,
        code_params,
        cppargs,
    )

    compiled_ode.field_params = {}
    compiled_ode.changed_scalar_parameters = {}
    compiled_ode.changed_field_parameters = []
    compiled_ode.set_parameter = types.MethodType(
        _set_parameter,
        compiled_ode,
        _gosscpp.ParameterizedODE,
    )

    # Store gotran ode model
    compiled_ode._gotran = ode
    compiled_ode._field_parameters = field_parameters
    compiled_ode._field_states = field_states

    return compiled_ode


# dolfin_jit.func_doc = goss_jit.__doc__


def family_and_degree_from_str(space):
    assert isinstance(space, str), "expected a str as the 'space' argument"

    space = space.split("_")
    assert all(space) and len(space) == 2, (
        "Expected a family name (CG, "
        "DG, P, Quadrature) and a degree separated by '_' as the "
        "space argument"
    )

    if space[0] in ["Lagrange", "CG", "P"]:
        assert space[1] == "1", "expected only P1 spaces"
        family = "Lagrange"
        degree = 1
    elif space[0] == "Quadrature":
        assert space[1] in [
            str(d) for d in range(1, 10)
        ], "expected only Quadrature spaces with degree 1-9."
        family = "Quadrature"
        degree = int(space[1])
    elif space[0] in ["DG", "Discontinuous Lagrange"]:
        assert space[1] in [str(d) for d in range(10)], (
            "expected only " "discontinuous spaces with degree 0-9."
        )
        family = "Discontinuous Lagrange"
        degree = int(space[1])
    else:
        assert False, (
            "Expected a family name (CG, "
            "DG, P, Quadrature) and a degree separated by '_' as the "
            "space argument"
        )

    return family, degree


class DOLFINODESystemSolver(object):
    """
    DOLFINODESystemSolver is an adapter class for goss.ODESystemSolver
    making interfacing DOLFIN easier
    """

    @staticmethod
    def default_parameters():

        # Include default params for ImplicitODESolvers
        # FIXME: Considre adding more?
        goss_solver_params = ImplicitODESolver.default_parameters()
        solver_params = d.Parameters("solver_params")

        for key, value in goss_solver_params.items():
            solver_params.add(key, value)

        params = d.Parameters("DOLFINODESystemSolver")
        params.add("solver", "ImplicitEuler", goss_solvers)
        params.add("num_threads", 0, 0, 100)
        params.add(solver_params)
        params.add("use_cuda", False)

        if enable_cuda:
            default_cuda_params = goss.cuda.CUDAODESystemSolver.default_parameters()
            cuda_params = d.Parameters("cuda_params")
            cuda_params.add("block_size", 256, 1, 1024)
            cuda_params.add("nvcc_options", "")
            cuda_params.add(
                "solver",
                default_cuda_params["solver"],
                dict.get(default_cuda_params, "solver")._options,
            )
            cuda_params.add("float_precision", "double", ["single", "double"])
            cuda_params.add("ldt", default_cuda_params["ldt"], -1.0, 1.0e6)
            params.add(cuda_params)

        return params

    def __init__(  # noqa: C901
        self,
        mesh,
        odes,
        domains=None,
        space="P_1",
        params=None,
    ):
        """
        Arguments:
        ----------
        mesh : Mesh
           The domain. By default ODE dofs lives on the vertices.
        odes : ParameterizedODE, dict
           The ParameterizedODE, to be solved in each mesh entity. If a dict the
           keys defines the domain and the values defines the ODE should be
           solved on.
        domains : MeshFunction (optional)
           A MeshFunction describing the distribution of the ODEs over the mesh.
        space : str (optional)
           A str representing the space the dolfin Function holding fields
           states should live in. For now only P1, DGX, and QX are allowed.
        params : dict (optional)
           A dict providing parameters for the Solvers
        """
        # FIXME: This function is way to long
        from . import _gosscpp

        assert isinstance(mesh, (d.Mesh, d.Domain)), (
            "expected a dolfin Mesh " "or domain for the mesh argument"
        )
        assert isinstance(odes, (dict, _gosscpp.ParameterizedODE)), (
            "expected a" " dict or a ParametersizedODE for the odes argument"
        )

        # Get family and degree from str
        family, degree = family_and_degree_from_str(space)

        params = params or {}

        self.parameters = self.default_parameters()
        self.parameters.update(params)

        top_dim = mesh.topology().dim()

        if isinstance(odes, _gosscpp.ParameterizedODE) or len(odes) == 1:
            # Get rid of any default labels.
            if isinstance(odes, dict):
                odes = odes.values()[0]

            num_field_states = odes.num_field_states()
            field_names = odes.get_field_state_names()
            odes = {0: odes}
            distinct_domains = [0]
            assert domains is None, (
                "domains only expected when more than 1 " "ODE is given"
            )
        else:
            assert isinstance(domains, d.MeshFunctionSizet), (
                "expected a "
                "MeshFunction as the domains argument when more than "
                "1 ODE is given"
            )
            expected_dim = 0 if family == "Lagrange" else top_dim
            assert (
                domains.dim() == expected_dim
            ), "expected a domain to be a " "MeshFunction of topological dimension {} for {} space".format(
                expected_dim,
                space,
            )

            # Check given domains
            distinct_domains = list(sorted(set(domains.array())))
            assert d.MPI.max(mesh.mpi_comm(), len(distinct_domains)) == d.MPI.max(
                mesh.mpi_comm(),
                len(odes),
            ), (
                "expected the number "
                "of distinct domains to be the same as the number of ODEs"
            )

            # Check and compare the number of field states
            ode_list = odes.values()
            last_ode = ode_list.pop()
            assert all(
                last_ode.num_field_states() == ode.num_field_states()
                for ode in ode_list
            ), ("expected all odes to have the " "same number of field states")

            last_field_state_names = last_ode.get_field_state_names()
            assert all(
                last_field_state_names == ode.get_field_state_names()
                for ode in ode_list
            ), (
                "expected all odes to have the same name and order "
                "of the field states (Might be changed in the future.)"
            )

            num_field_states = last_ode.num_field_states()
            field_names = last_ode.get_field_state_names()

        # Create dof storage for the dolfin function
        if num_field_states > 1:
            V = d.VectorFunctionSpace(mesh, family, degree, dim=num_field_states)
        else:
            V = d.FunctionSpace(mesh, family, degree)

        # Get the mesh entity to dof mappings
        mesh_entity_to_dof = entity_to_dofs(V)
        first_dof, last_dof = V.dofmap().ownership_range()
        num_local_dofs = last_dof - first_dof

        # Extract dof and mesh entity information
        dof_maps = {}
        goss_values = {}
        num_dofs = {}
        goss_indices = {}
        mesh_entities = {}

        nested_dofs = len(distinct_domains) > 1 or num_field_states > 1

        if enable_cuda and self.parameters.use_cuda:
            float_type = (
                np.float32
                if params["cuda_params"]["float_precision"] == "single"
                else np.float64
            )
        else:
            float_type = np.float_

        for label in distinct_domains:
            dof_maps[label] = OrderedDict((key, []) for key in field_names)
            if num_field_states > 1:
                num_entity_dofs_scalar = V.sub(0).dofmap().num_entity_dofs(top_dim)
            else:
                num_entity_dofs_scalar = V.dofmap().num_entity_dofs(top_dim)

            # If domains given
            if domains:
                local_mesh_entities = (domains.array() == label).nonzero()[0]
                if family == "Lagrange":
                    mesh_entities[label] = local_mesh_entities
                else:
                    mesh_entities_all = np.zeros(
                        len(local_mesh_entities) * num_entity_dofs_scalar,
                        dtype=np.uintp,
                    )
                    for index in range(num_entity_dofs_scalar):
                        mesh_entities_all[index::num_entity_dofs_scalar] = (
                            local_mesh_entities * num_entity_dofs_scalar + index
                        )

                    mesh_entities[label] = mesh_entities_all

            else:
                num_entities = (
                    mesh.num_vertices()
                    if family == "Lagrange"
                    else mesh.num_cells() * num_entity_dofs_scalar
                )
                mesh_entities[label] = np.arange(num_entities, dtype=np.uintp)

            for local_dof_offset in range(num_field_states):
                local_entities = (
                    mesh_entities[label] * num_field_states + local_dof_offset
                )
                dofs = mesh_entity_to_dof[local_entities]
                dof_maps[label][field_names[local_dof_offset]] = dofs[
                    (0 <= dofs) * (dofs < num_local_dofs)
                ]

            # Check we have the same number of dofs per field name.
            if len(field_names) > 1:
                assert all(
                    len(dof_maps[label][field_names[0]])
                    == len(dof_maps[label][field_name])
                    for field_name in field_names[1:]
                ), ("expected all " "fields to have the same number of dofs")

            # Num dofs per field state per label (only store one of the
            # field states as it has to be the same for all field states)
            num_dofs[label] = len(dof_maps[label][field_names[0]])

            # Store the dofs as numpy arrays
            for names, value in dof_maps[label].items():
                dof_maps[label][names] = np.array(value, dtype=np.intc)

            # Allocate memory for value transfer to ODESystemSolver
            goss_values[label] = np.zeros(
                num_dofs[label] * num_field_states,
                dtype=float_type,
            )

            # Create a nested index set for putting values into field_states
            # in a ODESystemSolver
            goss_indices[label] = OrderedDict(
                (
                    key,
                    np.arange(
                        offset,
                        num_dofs[label] * num_field_states + offset,
                        num_field_states,
                        dtype=np.intc,
                    ),
                )
                for offset, key in enumerate(field_names)
            )

        # Allocate memory for accessing values from DOLFIN Function
        if nested_dofs:

            # Always use np.float_ for the dolfin_values
            self._dolfin_values = np.concatenate(
                tuple(value for value in goss_values.values()),
            ).astype(np.float_)
        else:
            # If not nested_dofs just grab the goss_value array
            if float_type == np.float32:
                self._dolfin_values = goss_values.values()[0].astype(np.float_)
            else:
                self._dolfin_values = goss_values.values()[0]

        # A dof index array used to access the dofs from the dolfin vector
        dof_maps["dolfin"] = np.arange(len(self._dolfin_values), dtype=np.intc)

        # Store arguments
        self._mesh = mesh
        self._odes = odes
        self._domains = domains
        self._float_type = float_type
        self._state_space = V
        self._saved_states = {}

        if enable_cuda and self.parameters.use_cuda:

            default_cuda_params = goss.cuda.CUDAODESystemSolver.default_parameters()
            cuda_params = params["cuda_params"]

            # Transfer cuda_params to CUDAODESystemSolver parameters
            default_cuda_params["block_size"] = cuda_params["block_size"]
            default_cuda_params["nvcc_options"] = cuda_params["nvcc_options"].split()
            default_cuda_params["solver"] = cuda_params["solver"]
            default_cuda_params["ldt"] = cuda_params["ldt"]
            default_cuda_params["code"]["float_precision"] = cuda_params[
                "float_precision"
            ]

            cuda_parameters = {}
            for label, ode in odes.items():
                cuda_parameters[label] = default_cuda_params.copy()
                cuda_parameters[label]["code"]["parameters"][
                    "field_parameters"
                ] = ode._field_parameters or [""]
                cuda_parameters[label]["code"]["states"][
                    "field_states"
                ] = ode._field_states or [""]

            # Instantiate the ODESystemSolvers
            self._ode_system_solvers = OrderedDict(
                (
                    label,
                    goss.cuda.CUDAODESystemSolver(
                        num_dofs[label],
                        odes[label]._gotran,
                        params=cuda_parameters[label],
                    ),
                )
                for label in distinct_domains
            )

        else:

            # Instantiate the solver
            solver = eval(self.parameters["solver"], _gosscpp.__dict__, {})()

            for param, value in self.parameters["solver_params"].items():
                if param in solver.parameters:
                    solver.parameters[param] = value

            # Instantiate the ODESystemSolvers
            self._ode_system_solvers = OrderedDict(
                (
                    label,
                    _gosscpp.ODESystemSolver(
                        num_dofs[label],
                        solver.copy(),
                        odes[label],
                    ),
                )
                for label in distinct_domains
            )

            # Set num of threads
            for solver in self._ode_system_solvers.values():
                solver.set_num_threads(self.parameters["num_threads"])

        # Check for field parameters and set any changed parameters
        self._field_params_dofs = {}
        self._param_values = {}
        for label in distinct_domains:
            ode = odes[label]

            # Scalar parameters are set and we need to update the
            # CUDAODESystemSolver
            if (
                enable_cuda
                and self.parameters.use_cuda
                and ode.changed_scalar_parameters
            ):

                # Get CUDAODESystemSolver and set parameters
                cuda_solver = self._ode_system_solvers[label]
                cuda_solver.set_parameters(**ode.changed_scalar_parameters)

                # Reset changed scalar parameters
                ode.changed_scalar_parameters = {}

            # If there are no field parameters
            if ode.num_field_parameters() == 0:
                continue

            # Check for any field params on this domain
            if len(ode.field_params) > 0:

                # function space of parameter Function
                V_param = ode.field_params.values()[0].function_space()

                assert (
                    V_param.ufl_element().family() == family
                    and V_param.ufl_element().degree() == degree
                ), (
                    "Expected space of ODESystemSolver to be the same as "
                    "passed field_parameters."
                )

                # Mesh entity and dof mappings
                first_dof_param, last_dof_param = V_param.dofmap().ownership_range()
                num_local_param_dofs = last_dof_param - first_dof_param

                mesh_entity_to_dof_param = entity_to_dofs(V_param)

                # Get dofs local to geometric domain
                dofs = mesh_entity_to_dof_param[mesh_entities[label]]
                dofs = dofs[(0 <= dofs) * (dofs < num_local_param_dofs)]
                self._field_params_dofs[label] = dofs

            # Init memory for setting field_parameters
            self._param_values[label] = np.zeros(
                num_dofs[label] * ode.num_field_parameters(),
                dtype=float_type,
            )

            for local_id, param in enumerate(ode.get_field_parameter_names()):
                if param in ode.changed_field_parameters:
                    if nested_dofs:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters()
                        ] = (
                            ode.field_params[param]
                            .vector()
                            .array()[self._field_params_dofs[label]]
                        )
                    else:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters()
                        ] = (ode.field_params[param].vector().array())

                else:
                    self._param_values[label][
                        local_id :: ode.num_field_parameters()
                    ] = ode.get_parameter(param)

            # Reset any changed field parameters
            ode.changed_field_parameters = []

            # Set field_parameter values
            if label in self._param_values:
                self._ode_system_solvers[label].set_field_parameters(
                    self._param_values[label],
                )

        # Store dof mapping and field value storages
        self._dof_maps = dof_maps
        self._goss_values = goss_values
        self._goss_indices = goss_indices
        self._num_dofs = num_dofs
        self._field_names = field_names
        self._num_field_states = len(field_names)
        self._num_distinct_domains = len(distinct_domains)
        self._distinct_domains = distinct_domains
        self._nested_dofs = nested_dofs

    def update_parameters(self):
        """
        Update the values in the any changed parameters
        """
        for label in self._distinct_domains:
            ode = self._odes[label]
            # field_params_changed = False
            for local_id, param in enumerate(ode.get_field_parameter_names()):
                if param in ode.changed_field_parameters:
                    # field_params_changed = True
                    if self._nested_dofs:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters()
                        ] = (
                            ode.field_params[param]
                            .vector()
                            .array()[self._field_params_dofs[label]]
                        )
                    else:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters()
                        ] = (ode.field_params[param].vector().array())

            # Update system solver
            if ode.changed_field_parameters:
                self._ode_system_solvers[label].set_field_parameters(
                    self._param_values[label],
                )
                ode.changed_field_parameters = []

            # Scalar parameters are set and we need to update the
            # CUDAODESystemSolver
            if (
                enable_cuda
                and self.parameters.use_cuda
                and ode.changed_scalar_parameters
            ):

                # Get CUDAODESystemSolver and set parameters
                cuda_solver = self._ode_system_solvers[label]
                cuda_solver.set_parameters(**ode.changed_scalar_parameters)

                # Reset changed scalar parameters
                ode.changed_scalar_parameters = {}

    @property
    def field_names(self):
        return self._field_names

    @property
    def num_field_states(self):
        return len(self._field_names)

    @property
    def state_space(self):
        return self._state_space

    @property
    def num_distinct_domains(self):
        return self._num_distinct_domains

    def from_field_states(self, v):
        """
        Copy values in stored field states to v
        """
        assert isinstance(v, d.Function), "expected a Function as the 'v' argument"
        # FIXME: Add proper check!
        # assert v in self._state_space, "expected v to be in the state space"

        # Get values from dolfin
        values = self._dolfin_values

        # Update solver with new field_state values
        for label, ode_system_solver in self._ode_system_solvers.items():

            # Fetch solution from stored field states
            ode_system_solver.get_field_states(self._goss_values[label])

            # Iterate over the fields and collect values and put back
            # into dolfin transit array if nested dofs
            if self._nested_dofs:
                for field_name in self._field_names:
                    goss_indices = self._goss_indices[label][field_name]
                    dof_maps = self._dof_maps[label][field_name]

                    # Get each field for each distinct domain
                    values[dof_maps] = self._goss_values[label][goss_indices]

            elif self._float_type == np.float32:

                values[:] = self._goss_values[label]

        # Put solution back into DOLFIN Function
        v.vector()[self._dof_maps["dolfin"]] = values

    def to_field_states(self, v):
        """
        Copy values in v to stored field states
        """
        assert isinstance(v, d.Function), "expected a Function as the 'v' argument"
        # FIXME: Add proper check!
        # assert v in self._state_space, "expected v to be in the state space"

        # Get values from dolfin
        values = self._dolfin_values
        v.vector().get_local(values, self._dof_maps["dolfin"])

        # Update solver with new field_state values
        for label, ode_system_solver in self._ode_system_solvers.items():

            # Iterate over the fields and collect values if nested dofs
            if self._nested_dofs:
                for field_name in self._field_names:
                    goss_indices = self._goss_indices[label][field_name]
                    dof_maps = self._dof_maps[label][field_name]

                    # Get each field for each distinct domain
                    self._goss_values[label][goss_indices] = values[dof_maps]

            # If single precision we need to copy
            elif self._float_type == np.float32:

                self._goss_values[label][:] = values

            # Transfer values to Solver
            ode_system_solver.set_field_states(self._goss_values[label])

    def step(self, interval, v):
        """
        Solve on the given time step (t0, t1).

        End users are recommended to use solve instead.

        Arguments:
        interval : tuple
          The time interval (t0, t1) for the step
        """

        assert isinstance(v, d.Function), "expected a Function as the 'v' argument"
        # FIXME: Add proper check!
        # assert v in self._state_space, "expected v to be in the state space"

        timer = d.Timer("ODE step")  # noqa: F841
        (t0, t1) = interval
        dt = t1 - t0

        # Update local field states
        self.to_field_states(v)

        # Update any changed field_parameters
        self.update_parameters()

        # Step solvers
        for label, ode_system_solver in self._ode_system_solvers.items():
            ode_system_solver.forward(t0, dt)

        # Copy solution from local field states
        self.from_field_states(v)

    def save_states(self):
        """
        Save the present state
        """
        # If not create copy the states and exit
        if not self._saved_states:
            if enable_cuda and self.parameters.use_cuda:
                self._saved_states = dict(
                    (label, ode_system_solver.get_cuda_states().copy())
                    for label, ode_system_solver in self._ode_system_solvers.items()
                )

            else:
                self._saved_states = dict(
                    (label, ode_system_solver.states().copy())
                    for label, ode_system_solver in self._ode_system_solvers.items()
                )
            return

        # Save states for each solver
        if enable_cuda and self.parameters.use_cuda:
            for label, ode_system_solver in self._ode_system_solvers.items():
                self._saved_states[label][:] = ode_system_solver.get_cuda_states()

        else:
            for label, ode_system_solver in self._ode_system_solvers.items():
                self._saved_states[label][:] = ode_system_solver.states()

    def restore_states(self):
        """
        Restore the states from any saved states
        """
        # If not create copy the states and exit
        if not self._saved_states:
            raise RuntimeError("Cannot restore any states when none are saved.")

        # Restor the states from the saved one
        if enable_cuda and self.parameters.use_cuda:
            for label, ode_system_solver in self._ode_system_solvers.items():
                ode_system_solver.set_cuda_states(self._saved_states[label])

        else:
            for label, ode_system_solver in self._ode_system_solvers.items():
                ode_system_solver.states()[:] = self._saved_states[label]
