import os

import numpy as np
import pytest
from goss import goss_solvers
from goss import jit
from goss import ODESystemSolver
from goss.solvers import RL1
from gotran import load_ode


parametrize = pytest.mark.parametrize


def convergence_order(errors, base=2):
    import math

    orders = [0.0] * (len(errors) - 1)
    for i in range(len(errors) - 1):
        try:
            orders[i] = math.log(errors[i] / errors[i + 1], base)
        except ZeroDivisionError:
            orders[i] = np.nan

    return orders


class TestODESolvers(object):
    orders = dict(
        RK4=4,
        RK2=2,
        RL1=1,
        RL2=2,
        GRL1=1,
        GRL2=2,
        ExplicitEuler=1,
        ThetaSolver=2,
        BasicImplicitEuler=1,
        RKF32=3,
        ESDIRK23a=3,
        ESDIRK4O32=4,
        ImplicitEuler=1,
    )
    oscilator = jit(load_ode("oscilator"))
    exclude_osc = ["ESDIRK4O32", "ESDIRK23a"]
    exclude_tent = ["RKF32", "ESDIRK4O32", "ESDIRK23a", "BasicImplicitEuler"]
    dir_path = os.path.dirname(__file__)
    # Vm_reference = np.fromfile(os.path.join(dir_path, "Vm_reference.npy"))
    tentusscher = jit(load_ode("tentusscher_2004_mcell"))

    @parametrize(("solver_str"), goss_solvers)
    def test_convergence_order(self, solver_str):

        if solver_str in self.exclude_osc:
            return

        solver = eval(solver_str)(self.oscilator)

        tstop = 10.0
        exact = np.array([np.cos(tstop), np.sin(tstop)])
        errors = []
        for dt in [0.05, 0.025, 0.0125, 0.00625]:
            u = np.array([1.0, 0.0])
            t = 0.0
            nsteps = int(tstop / dt)
            for step in range(nsteps):
                solver.forward(u, t, dt)
                t += dt

            errors.append(np.sqrt(np.sum(((exact - u) / exact) ** 2)))

        assert min(convergence_order(errors)) >= self.orders[solver_str] - 0.1

    @parametrize(("solver_str"), goss_solvers)
    def test_long_run(self, solver_str):

        if solver_str in self.exclude_tent:
            return

        dt = 0.0002
        tstop = 10
        ind_V = self.tentusscher.get_state_names().index("V")
        # dt_ref = 0.1

        solver = eval(solver_str)(self.tentusscher)

        self.tentusscher.set_parameter("stim_amplitude", 52.0)
        self.tentusscher.set_parameter("stim_start", 0.0)

        u = self.tentusscher.init_state_values()

        t = 0.0
        nsteps = int(tstop / dt)
        for step in range(nsteps):
            solver.forward(u, t, dt)
            t += dt

        # Test against run with scipy integrate
        assert abs(u[ind_V] - 12.948) < 1e-3


class TestODESystemSolver(object):

    ode = jit(
        load_ode("tentusscher_2004_mcell"),
        field_states=["V", "Ca_i"],
        field_parameters=["g_CaL", "K_o"],
    )

    def test_ode_interface(self):
        assert self.ode.num_field_states() == 2
        assert self.ode.num_field_parameters() == 2
        assert self.ode.num_states() == 17
        assert self.ode.num_parameters() == 45

        Cm = self.ode.get_parameter("Cm")
        assert Cm == 0.185

        self.ode.set_parameter("Cm", 0.2)
        assert self.ode.get_parameter("Cm") == 0.2
        self.ode.set_parameter("Cm", 0.185)

        assert isinstance(self.ode.get_state_names(), list)
        assert isinstance(self.ode.get_parameter_names(), list)
        assert len(self.ode.get_state_names()) == 17
        assert len(self.ode.get_parameter_names()) == 45

        assert "stim_start" in self.ode.get_parameter_names()
        state_names = [
            "Xr1",
            "Xr2",
            "Xs",
            "m",
            "h",
            "j",
            "d",
            "f",
            "fCa",
            "s",
            "r",
            "g",
            "Ca_i",
            "Ca_SR",
            "Na_i",
            "V",
            "K_i",
        ]
        state_names.sort()
        ode_state_names = self.ode.get_state_names()
        ode_state_names.sort()
        assert ode_state_names == state_names

        with pytest.raises(RuntimeError):
            self.ode.set_parameter("cm", 0.185)

        with pytest.raises(RuntimeError):
            self.ode.get_parameter("JADA")

    def test_ode_system_interface(self):

        num_nodes = 100
        solver = RL1()
        system = ODESystemSolver(num_nodes, solver, self.ode)
        system.reset_default()
        field_states = np.zeros(self.ode.num_field_states() * num_nodes)
        # field_parameters = np.zeros(self.ode.num_field_parameters() * num_nodes)

        assert isinstance(system.states(), np.ndarray)
        assert len(system.states()) == num_nodes * self.ode.num_states()

        # Check tangled access
        system.get_field_states(field_states, True)
        assert sum(field_states[::2] == -86.2) == num_nodes
        assert sum(field_states[1::2] == 0.0002) == num_nodes

        # Check untangled access
        system.get_field_states(field_states, False)
        assert sum(field_states[:num_nodes] == -86.2) == num_nodes
        assert sum(field_states[num_nodes:] == 0.0002) == num_nodes

        assert system.num_nodes() == num_nodes

        with pytest.raises(ValueError):
            system.set_field_states(
                np.zeros(self.ode.num_field_states() * num_nodes // 2),
            )

        with pytest.raises(TypeError):
            system.set_field_states(
                np.zeros(self.ode.num_field_states() * num_nodes, dtype=int),
            )

        with pytest.raises(ValueError):
            system.set_field_parameters(
                np.zeros(self.ode.num_field_states() * num_nodes // 2),
            )

        with pytest.raises(TypeError):
            system.set_field_parameters(
                np.zeros(self.ode.num_field_states() * num_nodes, dtype=int),
            )

        # system.get_field_state_components()
