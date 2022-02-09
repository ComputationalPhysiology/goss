import abc


class ODESolver(abc.ABC):
    pass


class ExplicitEuler(ODESolver):
    pass


class RL1(ODESolver):
    pass


class ImplicitODESolver(ODESolver):
    pass
