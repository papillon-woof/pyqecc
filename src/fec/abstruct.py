from abc import ABCMeta, abstractmethod
import numpy as np

class CODE(metaclass=ABCMeta):
    _name = ""
    def __init__(self,n,k):
        self._n = n
        self._k = k
        self._nk = n - k
        self._R = self.k/self.n

    '''
    #未実装
    @property
    def enc_circ(self):
        return self._enc_circuit
    '''

    @property
    def n(self):
        return self._n

    @property
    def k(self):
        return self._k

    @property
    def R(self):
        return self._R

    @property
    def name(self):
        return self._name

    def __str__(self):
        output = ""
        output+="codename        :"+str(self.name)+"\n"
        output+="n               : "+str(self.n)+"\n"
        output+="k               : "+str(self.k)+"\n"
        output+="R               : "+str(self.R)+"\n"
        return output
