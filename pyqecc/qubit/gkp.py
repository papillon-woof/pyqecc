import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from .abstruct import *
# measure
Measure = None
SW_gate = None
# 1 qubit gate
# (p,q)のシフト量
M_gate = np.array([[1,0],[0,1]])
I_gate = np.array([[1,0],[0,1]])
X_gate = np.array([[1,0],[0,np.sqrt(np.pi)]])
Z_gate = np.array([[np.sqrt(np.pi),0],[0,1]])
Y_gate = np.array([[np.sqrt(np.pi),0],[[0,np.sqrt(np.pi)],0]])
H_gate = np.array([[0,1],[-1,0]])
S_gate = np.array([[1,-1],[0,1]])

# 2 qubit gate (p0,p1,...pn,q0,q1,...,qn)
CX_gate = np.array([[1,0,1,0],[0,1,0,0],[0,0,1,0],[0,-1,0,1]])
#CZ_gate = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]],dtype=np.complex128)

class GKP_QUBIT(Qubits):
    _name = "qubits"
    gate_fundamental_dictionary = {
        "M":M_gate,
        "X":X_gate,
        "Y":Y_gate,
        "Z":Z_gate,
        "T":T_gate,
        "CX":CX_gate,
        #"CZ":CZ_gate,
        "SW":SW_gate,
        "I":I_gate,
    }
    def __init__(self,n,delta=0.2,type="square",photon_number=5,valid=10,bit_bin=0.01):
        self.bit_bin = bit_bin
        self.delta = delta
        self._n = n
        self.length = 2**valid
        self.wigner = np.zeros((self.n,self.length*2+1))
        self.photon_number = photon_number

        # shift value for GKP qubits
        self.shift_value = np.np.zeros(2 * self.n)
        self.wigner_function()

    def SW(self,name,first_idx=-1,second_idx=-1):
        if first_idx<self.n:
            raise ValueError("n<idx")
        self.shift_value[2*first_idx],self.shift_value[2*second_idx] = self.shift_value[2*first_idx],self.shift_value[2*second_idx]
        self.shift_value[2*first_idx+1],self.shift_value[2*second_idx+1] = self.shift_value[2*first_idx],self.shift_value[2*second_idx+1]

    def gate(self,name,first_idx=-1,second_idx=-1):
        self.shift_value = np.dots(self.gate_fundamental_dictionary[name],self.shift_value)
        self.shift_value = np.dots(self.gate_fundamental_dictionary[name],self.shift_value)

    def wigner_function(self):
        q = self.bit_bin*np.array([i for i in range(-self.length,self.length+1)])
        for num in range(self.n):
            for s in range(-self.photon_number,self.photon_number+1):
                for j in range(-self.length,self.length+1):
                    self.wigner[num,j] += np.exp(-(self.delta**2)/2*((2*s)**2)*np.pi)*np.exp(-(1/(2*self.delta**2)*((q[j]-2*s*np.sqrt(np.pi))**2)))
        plt.plot(self.wigner[0])
        plt.show()

    @property
    def n(self):
        return self._n