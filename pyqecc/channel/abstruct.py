from abc import ABCMeta, abstractmethod
import numpy as np

class Channel(metaclass=ABCMeta):
    _name = ""

    def __init__(self,seed):
        if not seed is None:
            np.random.seed(seed=seed)
    
    @abstractmethod
    def channel(self):
        pass