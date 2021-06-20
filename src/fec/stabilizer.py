import numpy as np
from .abstruct import *
from ..util import *
class SC(CODE):
    def __init__(self,n,k,H='random',beta=None):
        self._name = "stabilizer"
        self._n = n
        self._k = k
        self._R = self._k/self._n
        if H in ['random']:
            pass
        else:
            self._H = H
        self._syndrome = beta

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
    def get_syndrome(self,e):
        return symplex_binary_inner_product(self._H,e)

    def hard_decode(self,e):
        s = symplex_binary_inner_product(self._H,e)
        ee = np.zeros(2*self.n,dtype='i1')
        for i in range(self.n-self.k):
            ee=s[i]*self._syndrome[i]+ee
        ee = np.mod(ee,2)
        return ee

    def __str__(self):
        output = ""
        output+="codename        :"+str(self.name)+"\n"
        output+="n               : "+str(self.n)+"\n"
        output+="k               : "+str(self.k)+"\n"
        output+="R               : "+str(self.R)+"\n"
        return output
