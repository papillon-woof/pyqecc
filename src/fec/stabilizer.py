import numpy as np
from .abstruct import *
class SC(CODE):
    def __init__(self,n,k,Htype='random'):
        self._n = n
        self._k = k
        self._H = None
        self.enc_circuit = None
        self.dec_circuit = None

    #量子情報ビット
    def get_enc_circ(self):
        return self.enc_circuit

    #シンドロームと確率分布を受け取る．
    def get_dec_circ(self):
        return self.dec_circuit

    # soft decision
    def ml_dec(self,p):
        return L

    # soft decision
    def rn_dec(self,p):
        return L

    # hard decision
    def syndrome_dec(self,e):
        return symplex_binary_inner_product(H,e)
        return E
