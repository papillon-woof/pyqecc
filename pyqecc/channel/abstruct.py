import itertools
from abc import ABCMeta, abstractmethod
import numpy as np
class Channel(metaclass=ABCMeta):
    _name = ""

    def __init__(self,seed,n=-1):
        self._channel_parameter_name = [] # string
        self._channel_parameter = {} # list or np.array
        self._channel_output = {"E",None}
        self._n = n
        if not seed is None:
            np.random.seed(seed=seed)

    def generate_param(self):
        p = []
        for v in self._channel_parameter.values():
            p = itertools.product(v, p)
        self._channel_parameter_name = self._channel_parameter.keys()
        self._channel_parameter["PARAM_SET"] = list(set(p))


    @abstractmethod
    def channel(self):
        return self.channel_output
        #Need a "return self.channel_output"
    
    def set_n(self,n):
        self._n = n
        self._channel_output["E"] = np.zeros(n,dtype='i1')
    
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