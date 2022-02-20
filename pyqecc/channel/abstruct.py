import itertools
from abc import ABCMeta, abstractmethod
import numpy as np


class Channel(metaclass=ABCMeta):
    _name = ""

    def __init__(self, n, seed):
        self._channel_parameter_name = []  # string
        self._channel_parameter = {}  # list or np.array
        self._channel_output = {"E": None}
        self.n = n
        if not seed is None:
            np.random.seed(seed=seed)

    def generate_param(self):
        key = list(self._channel_parameter.keys())
        p = self._channel_parameter[key[0]]
        for k in key[1:]:
            p = itertools.product(self._channel_parameter[k], p)
        p = [v for v in p]
        self._channel_parameter_name = self._channel_parameter.keys()
        self._channel_parameter["PARAM_SET"] = p

    @abstractmethod
    def channel(self):
        return self.channel_output
        # Need a "return self.channel_output"

    def set_n(self, n):
        self._n = n
        self._channel_output["E"] = np.zeros(2 * n, dtype="i1")

    @property
    def channel_parameter_name(self):
        return self._channel_parameter_name

    @property
    def channel_parameter(self):
        return self._channel_parameter

    @property
    def channel_output(self):
        return self._channel_output

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, n):
        if n <= 0:
            raise ValueError("PHYSICAL_QUBIT is more than 0")
        self._n = n
        self._channel_output["E"] = np.zeros(2 * n, dtype="i1")
