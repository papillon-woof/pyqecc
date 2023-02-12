from typing import List, Union, Tuple
from abc import ABCMeta, abstractmethod
import numpy as np
import numpy.typing as npt


class CODE(metaclass=ABCMeta):
    _name = ""

    def __init__(self, n:int, k:int):
        self._n = n
        self._k = k
        self._nk = n - k
        self._R = self.k / self.n

    """
    #not implimented
    @property
    def enc_circ(self):
        return self._enc_circuit
    """

    @abstractmethod
    def get_syndrome(self) -> Union[float,List[float],npt.NDArray[np.complex64],Tuple[float,float]]:
        pass

    @abstractmethod
    def decode(self):
        pass

    @property
    def n(self):
        return self._n

    @property
    def k(self):
        return self._k

    @property
    def nk(self):
        return self._nk

    @property
    def R(self):
        return self._R

    @property
    def name(self):
        return self._name

    def __str__(self):
        output = ""
        output += "CODE_NAME       :" + str(self.name) + "\n"
        output += "PHYSICAL_QUBITS : " + str(self.n) + "\n"
        output += "LOGICAL_QUBITS  : " + str(self.k) + "\n"
        output += "CODE_RATE       : " + str(self.R) + "\n"
        return output
