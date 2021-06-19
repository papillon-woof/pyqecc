import numpy as np

class SC(CODE):
    def __init__(self,n,k,Htype='random'):
        self._n = n
        self._k = k
        self._H = None

    #量子情報ビット
    def enc(self):
        return quantum_cercuit

    #シンドロームと確率分布を受け取る．
    def dec(self):
        return quantum_cercuit

    # soft decision
    def ml_dec(self,p):
        return L

    # soft decision
    def rn_dec(self,p):
        return L

    # hard decision
    def syndrome_dec(self,p):
        return E
