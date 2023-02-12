from typing import Dict, List
import itertools
from abc import ABCMeta, abstractmethod
import numpy as np


class Channel(metaclass=ABCMeta):
    _name = ""

    def __init__(self, n: int, seed: int) -> None:
        self._channel_parameter_name: List[str] = []  # string
        self._channel_parameter: Dict = {}  # list or np.array
        self._channel_output: Dict = {"E": None}
        self.n = n
        if seed != -1:
            np.random.seed(seed=seed)

    def generate_param(self) -> None:
        key = list(self._channel_parameter.keys())
        p = self._channel_parameter[key[0]]
        for k in key[1:]:
            p = itertools.product(self._channel_parameter[k], p)
        p = [v for v in p]
        self._channel_parameter_name = list(self._channel_parameter.keys())
        self._channel_parameter["PARAM_SET"] = p

    @abstractmethod
    def channel(self) -> Dict:
        return self.channel_output
        # Need a "return self.channel_output"

    def set_n(self, n: int) -> None:
        self._n = n
        self._channel_output["E"] = np.zeros(2 * n, dtype="i1")

    @property
    def channel_parameter_name(self) -> List:
        return self._channel_parameter_name

    @property
    def channel_parameter(self) -> Dict:
        return self._channel_parameter

    @property
    def channel_output(self) -> Dict:
        return self._channel_output

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, n: int):
        if n <= 0:
            raise ValueError("PHYSICAL_QUBIT is more than 0")
        self._n = n
        self._channel_output["E"] = np.zeros(2 * n, dtype="i1")
