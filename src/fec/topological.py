import numpy as np
from .abstruct import *
from ..util import *

class TORIC(SC):
    def __init__(self,n,k,G='random',T=None,L=None,P=None,iid=True):
        #H is replesentation for graph into adjacent matrix:
        self._G = G
        self._H = g2c(H)
        self._name = "toric code"
        self._n = n
        self._k = k
        self._R = self._k/self._n
        if H in ['random']:
            pass
        else:
            self._H = H
        self._T = T
        self._L = L

        self.enc_circuit = None
        self.dec_circuit = None
        self._P = self.set_P(P)
        self.ML_decoding_qubit_limit = 15
        
    def decoder(self,H):
        pass
