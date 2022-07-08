import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from .abstruct import *

class GKP_QUBIT(Qubits):
    def __init__(self,n,delta=0.2,type="square",photon_number=5,valid=10,bit_bin=0.01):
        self.bit_bin = bit_bin
        self.delta = delta
        self._n = n
        self.length = 2**valid
        self.wigner = np.zeros((self.n,self.length*2+1))
        self.photon_number = photon_number
        
        self.wigner_function()

    def wigner_function(self):
        q = self.bit_bin*np.array([i for i in range(-self.length,self.length+1)])
        for num in range(self.n):
            for s in range(-self.photon_number,self.photon_number+1):
                for j in range(-self.length,self.length+1):
                    self.wigner[num,j] += np.exp(-(self.delta**2)/2*((2*s)**2)*np.pi)*np.exp(-(1/(2*self.delta**2)*((q[j]-2*s*np.sqrt(np.pi))**2)))
        plt.plot(self.wigner[0])
        plt.show()

    def X(i):
        pass
        #ここに，シフト処理を書く
    
    def M(i):
        pass

    def T(i):
        pass
        #ここに，Tゲート処理を書く．

    @property
    def n(self):
        return self._n