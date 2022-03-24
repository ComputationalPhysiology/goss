from __future__ import annotations

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

__all__ = ["DOLFINODESystemSolver", "family_and_degree_from_str"]

import numpy as np

import typing

try:
    import dolfin
except ImportError as e:
    raise ImportError("dolfin is not present") from e

# Import Gotran and try import cuda solver
# from goss import goss_solvers
# from goss.solvers import ImplicitODESolver
from .ode import ParameterizedODE
from . import solvers
from .systemsolver import ODESystemSolver

# import goss.cuda

# enable_cuda = goss.cuda.cuda is not None
enable_cuda = False

KeyType = typing.Union[str, int]


def entity_to_dofs(V):
    assert isinstance(V, dolfin.FunctionSpace)
    mesh = V.mesh()
    dim = mesh.topology().dim()
    dm = V.dofmap()
    if V.ufl_element().family() == "Lagrange":
        return dolfin.vertex_to_dof_map(V)

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

    for cell in dolfin.cells(mesh):
        index = cell.index()
        for field_ind, dm_l in enumerate(dms):
            entity_to_dof[
                index * num_entity_dofs
                + field_ind : (index + 1) * num_entity_dofs
                + field_ind : num_fields
            ] = dms[field_ind].cell_dofs(index)

    return entity_to_dof


class DOLFINParameterizedODE(ParameterizedODE):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.field_params: dict[KeyType, dolfin.Function] = kwargs.get(
            "field_params",
            {},
        )
        self.changed_scalar_parameters: dict[KeyType, float] = kwargs.get(
            "changed_scalar_parameters",
            {},
        )
        self.changed_field_parameters: list[KeyType] = kwargs.get(
            "changed_field_parameters",
            [],
        )
        self._initial_conditions = self.default_initial_conditions()
        init_conditions = kwargs.get("init_conditions", {})
        self.set_initial_conditions(**init_conditions)

    def default_initial_conditions(self) -> dict[str, float]:
        return dict(zip(self.state_names, self.get_ic()))

    @classmethod
    def from_ode(cls, ode: ParameterizedODE):
        if isinstance(ode, DOLFINParameterizedODE):
            return ode
        field_params = getattr(ode, "field_params", {})
        changed_scalar_parameters = getattr(ode, "changed_scalar_parameters", {})
        changed_field_parameters = getattr(ode, "changed_field_parameters", [])
        return cls(
            ode._cpp_object,
            field_params=field_params,
            changed_scalar_parameters=changed_scalar_parameters,
            changed_field_parameters=changed_field_parameters,
        )

    def num_states(self):
        # Make compatible with cbcbeat
        return super().num_states

    def set_parameter(self, name, value):
        assert isinstance(name, str), "expected a str as the name argument"
        assert isinstance(
            value,
            (int, float, dolfin.Function),
        ), "expected a scalar or a Function for the value argument"
        if isinstance(value, (float, int)):
            self._cpp_object.set_parameter(name, value)
            self.changed_scalar_parameters[name] = value

        else:
            field_param_names = self.field_parameter_names
            assert (
                name in field_param_names
            ), f"'{name}' is not a field parameter in this ode"
            self.field_params[name] = value
            self.changed_field_parameters.append(name)

    def convert(
        odes: typing.Union[ParameterizedODE, dict[int, ParameterizedODE]],
    ) -> dict[int, DOLFINParameterizedODE]:
        if isinstance(odes, ParameterizedODE):
            return {0: DOLFINParameterizedODE.from_ode(odes)}

        return {
            label: DOLFINParameterizedODE.from_ode(ode) for label, ode in odes.items()
        }

    def set_initial_conditions(self, **init):
        "Update initial_conditions in model"
        for init_name, init_value in init.items():
            if init_name not in self._initial_conditions:
                raise RuntimeError("'{init_name}' is not a parameter in {self}")
            if not isinstance(init_value, (float, int)) and not isinstance(
                init_value._cpp_object,
                dolfin.cpp.function.GenericFunction,
            ):
                raise RuntimeError("'{init_name}' is not a scalar or a GenericFunction")
            if (
                hasattr(init_value, "_cpp_object")
                and isinstance(
                    init_value._cpp_object,
                    dolfin.cpp.function.GenericFunction,
                )
                and init_value._cpp_object.value_size() != 1
            ):
                raise RuntimeError("expected the value_size of '{init_name}' to be 1")
            self._initial_conditions[init_name] = init_value

    def initial_conditions(self) -> dolfin.Expression:
        "Return initial conditions for v and s as an Expression."
        return dolfin.Expression(
            list(self.state_names), degree=1, **self._initial_conditions
        )

    def field_state_initial_conditions(self) -> dolfin.Expression:
        "Return initial conditions for v and s as an Expression."
        ic = {
            name: value
            for name, value in self._initial_conditions.items()
            if name in self.field_state_names
        }
        if len(ic) == 1:
            # Make sure expression has rank 0
            return dolfin.Expression(*self.field_state_names, degree=1, **ic)
        # Otherwise expression will have rank 1
        return dolfin.Expression(self.field_state_names, degree=1, **ic)


def family_and_degree_from_str(space: str) -> typing.Tuple[str, int]:
    assert isinstance(space, str), "expected a str as the 'space' argument"

    family, degree = space.split("_")

    if family in ["P", "CG"]:
        family = "Lagrange"
    return family, int(degree)


class GossDofs(typing.NamedTuple):
    num_dofs: dict[KeyType, int]
    goss_values: dict[KeyType, np.ndarray]
    goss_indices: dict[KeyType, dict[KeyType, np.ndarray]]
    dof_maps: dict[KeyType, np.ndarray]
    mesh_entities: dict[KeyType, np.ndarray]
    nested_dofs: bool
    dolfin_values: np.ndarray


def first_value(d: dict):
    """Get the first value in a dictionary"""
    return next(iter(d.values()))


def setup_dofs(
    V: dolfin.FunctionSpace,
    field_names: list[str],
    domains: typing.Optional[dolfin.cpp.mesh.MeshFunctionSizet] = None,
    distinct_domains: typing.Optional[list[int]] = None,
):

    mesh = V.mesh()

    num_field_states = len(field_names)
    # TODO assert size of V matches num_field_states

    family = V.ufl_element().family()

    top_dim = mesh.topology().dim()
    if distinct_domains is None:
        distinct_domains = [0]
        if domains is not None:
            check_domains_dim(domains, family)
            distinct_domains = list(sorted(set(domains.array())))
    # Create dof storage for the dolfin function

    # Get the mesh entity to dof mappings
    mesh_entity_to_dof = entity_to_dofs(V)
    first_dof, last_dof = V.dofmap().ownership_range()
    num_local_dofs = last_dof - first_dof

    # Extract dof and mesh entity information
    num_dofs: dict[KeyType, int] = {}
    goss_values: dict[KeyType, np.ndarray] = {}
    goss_indices: dict[KeyType, dict[KeyType, np.ndarray]] = {}
    dof_maps: dict[KeyType, np.ndarray] = {}
    mesh_entities: dict[KeyType, np.ndarray] = {}

    nested_dofs = len(distinct_domains) > 1 or num_field_states > 1

    float_type = np.float_

    for label in distinct_domains:
        dof_maps[label] = {key: [] for key in field_names}
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
            local_entities = mesh_entities[label] * num_field_states + local_dof_offset
            dofs = mesh_entity_to_dof[local_entities]
            dof_maps[label][field_names[local_dof_offset]] = dofs[
                (0 <= dofs) * (dofs < num_local_dofs)
            ]

        # Check we have the same number of dofs per field name.
        if len(field_names) > 1:
            assert all(
                len(dof_maps[label][field_names[0]]) == len(dof_maps[label][field_name])
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
        goss_indices[label] = {
            key: np.arange(
                offset,
                num_dofs[label] * num_field_states + offset,
                num_field_states,
                dtype=np.intc,
            )
            for offset, key in enumerate(field_names)
        }

    # Allocate memory for accessing values from DOLFIN Function
    if nested_dofs:

        # Always use np.float_ for the dolfin_values
        dolfin_values = np.concatenate(
            tuple(value for value in goss_values.values()),
        ).astype(np.float_)
    else:
        # If not nested_dofs just grab the goss_value array
        if float_type == np.float32:
            dolfin_values = first_value(goss_values).astype(np.float_)
        else:
            dolfin_values = first_value(goss_values)

    # A dof index array used to access the dofs from the dolfin vector
    dof_maps["dolfin"] = np.arange(len(dolfin_values), dtype=np.intc)

    return GossDofs(
        num_dofs=num_dofs,
        goss_values=goss_values,
        goss_indices=goss_indices,
        dof_maps=dof_maps,
        mesh_entities=mesh_entities,
        nested_dofs=nested_dofs,
        dolfin_values=dolfin_values,
    )


def check_domains_dim(domains, family):
    top_dim = top_dim = domains.mesh().topology().dim()
    expected_dim = 0 if family == "Lagrange" else top_dim
    assert domains.dim() == expected_dim, (
        "expected a domain to be a "
        f"MeshFunction of topological dimension {expected_dim} "
        f"for space with family {family}, got {domains.dim()}"
    )


def check_domains(domains, odes, mesh, space) -> list[int]:

    family, degree = family_and_degree_from_str(space)

    assert isinstance(domains, dolfin.cpp.mesh.MeshFunctionSizet), (
        "expected a "
        "MeshFunction as the domains argument when more than "
        "1 ODE is given"
    )
    check_domains_dim(domains, family)

    # Check given domains
    distinct_domains = list(sorted(set(domains.array())))
    assert dolfin.MPI.max(mesh.mpi_comm(), len(distinct_domains)) == dolfin.MPI.max(
        mesh.mpi_comm(),
        len(odes),
    ), (
        "expected the number "
        "of distinct domains to be the same as the number of ODEs"
    )

    # Check and compare the number of field states
    first_ode = first_value(odes)
    assert all(
        first_ode.num_field_states == ode.num_field_states for ode in odes.values()
    ), "expected all odes to have the same number of field states"

    assert all(
        first_ode.field_state_names == ode.field_state_names for ode in odes.values()
    ), (
        "expected all odes to have the same name and order "
        "of the field states (Might be changed in the future.)"
    )
    return distinct_domains


class DOLFINODESystemSolver:
    """
    DOLFINODESystemSolver is an adapter class for goss.ODESystemSolver
    making interfacing DOLFIN easier
    """

    @staticmethod
    def default_parameters():
        return {"solver": "GRL1", "num_threads": 0, "use_cuda": False}

    @staticmethod
    def default_parameters_dolfin():
        p = dolfin.Parameters("DOLFINODESystemSolver")
        for k, v in DOLFINODESystemSolver.default_parameters().items():
            p.add(k, v)
        return p

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

        assert isinstance(mesh, dolfin.Mesh), (
            "expected a dolfin Mesh " "or domain for the mesh argument"
        )
        assert isinstance(odes, (dict, ParameterizedODE)), (
            "expected a" " dict or a ParametersizedODE for the odes argument"
        )

        # Get family and degree from str
        family, degree = family_and_degree_from_str(space)

        params = params or {}

        self.parameters = self.default_parameters()
        self.parameters.update(params)

        solver: solvers.ODESolver = eval(self.parameters["solver"], solvers.__dict__)()
        solver.update_parameters(self.parameters.get("solver_parameters", {}))

        odes = DOLFINParameterizedODE.convert(odes)

        distinct_domains = [0]
        if len(odes) > 1:
            distinct_domains = check_domains(domains, odes, mesh, space)
        # else:

        # FIXME: Do we really need to check this?
        # assert domains is None, (
        #     "domains only expected when more than 1 " "ODE is given"
        # )

        ode = first_value(odes)
        num_field_states = ode.num_field_states
        field_names = ode.field_state_names
        distinct_domains = list(odes.keys())

        if num_field_states > 1:
            V = dolfin.VectorFunctionSpace(mesh, family, degree, dim=num_field_states)
        else:
            V = dolfin.FunctionSpace(mesh, family, degree)

        dofs = setup_dofs(
            V=V,
            field_names=field_names,
            domains=domains,
            distinct_domains=distinct_domains,
        )

        float_type = np.float_
        # Store arguments
        self._mesh = mesh
        self._odes = odes
        self._domains = domains
        self._float_type = float_type
        self._state_space = V
        self._saved_states = {}

        # Current and previous solution
        self.vs = dolfin.Function(V)
        self.vs_ = dolfin.Function(V)

        # Instantiate the ODESystemSolvers
        self._ode_system_solvers = {
            label: ODESystemSolver(
                num_nodes=dofs.num_dofs[label],
                solver=solver.copy(),
                ode=odes[label],
            )
            for label in distinct_domains
        }
        # Set num of threads
        for solver in self._ode_system_solvers.values():
            solver.num_threads = self.parameters["num_threads"]

        # Check for field parameters and set any changed parameters
        self._field_params_dofs = {}
        self._param_values = {}
        for label in distinct_domains:
            ode = odes[label]

            # If there are no field parameters
            if ode.num_field_parameters == 0:
                continue

            # Check for any field params on this domain
            if len(ode.field_params) > 0:

                # function space of parameter Function
                V_param = first_value(ode.field_params).function_space()

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
                _dofs = mesh_entity_to_dof_param[dofs.mesh_entities[label]]
                _dofs = _dofs[(0 <= _dofs) * (_dofs < num_local_param_dofs)]
                self._field_params_dofs[label] = _dofs

            # Init memory for setting field_parameters
            self._param_values[label] = np.zeros(
                dofs.num_dofs[label] * ode.num_field_parameters,
                dtype=float_type,
            )

            for local_id, param in enumerate(ode.field_parameter_names):
                if param in ode.changed_field_parameters:
                    if dofs.nested_dofs:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters
                        ] = (
                            ode.field_params[param]
                            .vector()
                            .get_local()[self._field_params_dofs[label]]
                        )
                    else:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters
                        ] = (ode.field_params[param].vector().get_local())

                else:
                    self._param_values[label][
                        local_id :: ode.num_field_parameters
                    ] = ode.get_parameter(param)

            # Reset any changed field parameters
            ode.changed_field_parameters = []

            # Set field_parameter values
            if label in self._param_values:
                self._ode_system_solvers[label].field_parameters = (
                    self._param_values[label],
                )

        # Store dof mapping and field value storages
        self._dofs = dofs
        # self._dof_maps = dof_maps
        # self._goss_values = goss_values
        # self._goss_indices = goss_indices
        # self._num_dofs = num_dofs
        self._field_names = field_names
        self._num_field_states = len(field_names)
        self._num_distinct_domains = len(distinct_domains)
        self._distinct_domains = distinct_domains
        # self._nested_dofs = nested_dofs

        self.initial_conditions_to_field_states()
        self.from_field_states()
        self.vs_.assign(self.vs)

    def solution_fields(self):
        """
        Return tuple of previous and current solution objects.
        """
        return (self.vs_, self.vs)

    def update_parameters(self):
        """
        Update the values in the any changed parameters
        """
        for label in self._distinct_domains:
            ode = self._odes[label]
            # field_params_changed = False
            for local_id, param in enumerate(ode.field_parameter_names):
                if param in ode.changed_field_parameters:
                    # field_params_changed = True
                    if self._dofs.nested_dofs:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters
                        ] = (
                            ode.field_params[param]
                            .vector()
                            .get_local()[self._field_params_dofs[label]]
                        )
                    else:
                        self._param_values[label][
                            local_id :: ode.num_field_parameters
                        ] = (ode.field_params[param].vector().get_local())

            # Update system solver
            if ode.changed_field_parameters:
                self._ode_system_solvers[label].field_parameters = (
                    self._param_values[label],
                )

                ode.changed_field_parameters = []

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

    def from_field_states(self):
        """
        Copy values in stored field states to v
        """
        # Get values from dolfin
        values = self._dofs.dolfin_values

        # Update solver with new field_state values
        for label, ode_system_solver in self._ode_system_solvers.items():

            # Fetch solution from stored field states
            self._dofs.goss_values[label] = ode_system_solver.field_states.ravel()

            # Iterate over the fields and collect values and put back
            # into dolfin transit array if nested dofs
            if self._dofs.nested_dofs:
                for field_name in self._field_names:
                    goss_indices = self._dofs.goss_indices[label][field_name]
                    dof_maps = self._dofs.dof_maps[label][field_name]

                    # Get each field for each distinct domain
                    values[dof_maps] = self._dofs.goss_values[label][goss_indices]

            else:
                values[:] = self._dofs.goss_values[label]

        # Put solution back into DOLFIN Function
        self.vs.vector()[self._dofs.dof_maps["dolfin"]] = values

    def initial_conditions_to_field_states(self):
        """
        Copy values from initial conditions to stored field states
        """
        V = self.vs.function_space()
        # Update solver with new field_state values
        for label, ode_system_solver in self._ode_system_solvers.items():

            ic = ode_system_solver.ode.field_state_initial_conditions()
            values = dolfin.interpolate(ic, V).vector().get_local()

            # Iterate over the fields and collect values if nested dofs
            if self._dofs.nested_dofs:
                for field_name in self._field_names:
                    goss_indices = self._dofs.goss_indices[label][field_name]
                    dof_maps = self._dofs.dof_maps[label][field_name]

                    # Get each field for each distinct domain
                    self._dofs.goss_values[label][goss_indices] = values[dof_maps]

            # If single precision we need to copy
            else:
                self._dofs.goss_values[label][:] = values

            ode_system_solver.field_states = self._dofs.goss_values[label]

    def to_field_states(self):
        """
        Copy values in v to stored field states
        """

        # Get values from dolfin
        values = self._dofs.dolfin_values
        values[:] = self.vs.vector()[self._dofs.dof_maps["dolfin"]]

        # Update solver with new field_state values
        for label, ode_system_solver in self._ode_system_solvers.items():

            # Iterate over the fields and collect values if nested dofs
            if self._dofs.nested_dofs:
                for field_name in self._field_names:
                    goss_indices = self._dofs.goss_indices[label][field_name]
                    dof_maps = self._dofs.dof_maps[label][field_name]

                    # Get each field for each distinct domain
                    self._dofs.goss_values[label][goss_indices] = values[dof_maps]

            # If single precision we need to copy
            else:
                self._dofs.goss_values[label][:] = values

            # Transfer values to Solver
            ode_system_solver.field_states = self._dofs.goss_values[label]

    def step(self, interval: typing.Tuple[float, float]) -> None:
        """Solve on the given time step (t0, t1).
        End users are recommended to use solve instead.

        Parameters
        ----------
        interval : typing.Tuple[float, float]
            The time interval (t0, t1) for the step
        v : dolfin.Function
            The function to store the solution
        """

        timer = dolfin.Timer("ODE step")  # noqa: F841
        (t0, t1) = interval
        dt = t1 - t0

        # Set-up current variables
        self.vs.assign(self.vs_)  # Start with good guess

        # Update local field states
        self.to_field_states()

        # Update any changed field_parameters
        self.update_parameters()

        # Step solvers
        for label, ode_system_solver in self._ode_system_solvers.items():
            ode_system_solver.forward(t0, dt)

        # Copy solution from local field states
        self.from_field_states()

    def save_states(self):
        """
        Save the present state
        """
        # If not create copy the states and exit
        if not self._saved_states:
            self._saved_states = dict(
                (label, ode_system_solver.states().copy())
                for label, ode_system_solver in self._ode_system_solvers.items()
            )
            return

        # Save states for each solver
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
        for label, ode_system_solver in self._ode_system_solvers.items():
            ode_system_solver.states()[:] = self._saved_states[label]
