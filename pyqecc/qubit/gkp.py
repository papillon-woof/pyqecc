import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from .abstruct import *
from ..channel.channel import GaussianQuantumChannel
# measure
Measure = None
SW_gate = None

# 1 qubit gate
# Shift (add or transformation)
M_gate = np.array([[1,0],[0,1]])
I_gate = np.array([[1,0],[0,1]])
X_gate = np.array([0,np.sqrt(np.pi)])
Z_gate = np.array([np.sqrt(np.pi),0])
Y_gate = np.array([[np.sqrt(np.pi),0],[[0,np.sqrt(np.pi)],0]])
H_gate = np.array([[0,1],[-1,0]])
S_gate = np.array([[1,-1],[0,1]])

# 2 qubit gate (p0,p1,...pn,q0,q1,...,qn)
CX_gate = np.array([[1,0,1,0],[0,1,0,0],[0,0,1,0],[0,1,0,-1]])

class GKPQubits(Qubits):
    _name = "qubits"
    gate_fundamental_dictionary = {
        "M":M_gate,
        "X":X_gate,
        "Y":Y_gate,
        "Z":Z_gate,
        "T":T_gate,
        "H":H_gate,
        "CX":CX_gate,
        #"CZ":CZ_gate,
        "SW":SW_gate,
        "I":I_gate,
    }
    def __init__(self,n,delta=0.2,type="square"):
        
        #################################
        # shift value for GKP qubits
        #################################
        self.delta = delta #finite squeezing
        self._n = n # qubits length
        self.shift_value = np.zeros(2 * self.n) #Shift value
        self.gaussian_quantum_channel = GaussianQuantumChannel(self.n,self.delta/np.sqrt(2)) # internal error model
        #################################              


    def SW(self,first_idx=-1,second_idx=-1,third_idx=-1,fourth_idx=-1):
        '''
        switching the 2 qubits.
        '''        
        if first_idx>self.n:
            raise ValueError("n<idx")
        if third_idx != -1 and fourth_idx != -1:
            # 1,2,3,4 <-> 2,1,4,3
            self.shift_value[2*first_idx],self.shift_value[2*second_idx],self.shift_value[2*third_idx],self.shift_value[2*fourth_idx] = self.shift_value[2*second_idx],self.shift_value[2*first_idx],self.shift_value[2*fourth_idx],self.shift_value[2*third_idx]
            self.shift_value[2*first_idx+1],self.shift_value[2*second_idx+1],self.shift_value[2*third_idx+1],self.shift_value[2*fourth_idx+1] = self.shift_value[2*second_idx+1],self.shift_value[2*first_idx+1],self.shift_value[2*fourth_idx+1],self.shift_value[2*third_idx+1]
        else:
            # 1,2 <-> 2,1
            self.shift_value[2*first_idx],self.shift_value[2*second_idx] = self.shift_value[2*second_idx],self.shift_value[2*first_idx]
            self.shift_value[2*first_idx+1],self.shift_value[2*second_idx+1] = self.shift_value[2*second_idx+1],self.shift_value[2*first_idx+1]
    
    def M(self,first_idx=-1,typeM="C"):#"C": get a continuos information "D": get a Z mesurement
        if "C" == typeM:
            return self.shift_value[2*first_idx:2*first_idx+2],
        elif "D" == typeM:
            self.SW(first_idx=0,second_idx=first_idx)
    
    def gate(self,name,first_idx=-1,second_idx=-1):
        if len(name)==2:
            self.SW(first_idx=0,second_idx=1,third_idx=first_idx,fourth_idx=second_idx)
            self.shift_value[:4] = np.dot(self.gate_fundamental_dictionary[name],self.shift_value[:4])
            self.SW(first_idx=0,second_idx=1,third_idx=first_idx,fourth_idx=second_idx)
        elif len(name)==1:
            if name == "M":
                return self.M(first_idx=first_idx)
            elif name == "X" or name == "Z":
                self.shift_value[2*first_idx:2*first_idx+2] += self.gate_fundamental_dictionary[name]

    def add_noise(self,u,v):
        if len(u)+len(v) == self.n * 2:
            ValueError("noise length is invalid.")
        err = np.c_[u,v].T.reshape(2 * self.n)
        self.shift_value = self.shift_value + err

    def internal_error(self):
        delta = self.gaussian_quantum_channel.channel(self.n)["shift_error"]
        self.add_noise(delta[:self.n],delta[self.n:])

    def get_shift_val(self,idx):
        return self.shift_value[2*idx],self.shift_value[2*idx+1]

    def wigner_function(self,photon_number=5,valid=10,bit_bin=0.01):
        length = 2 ** valid
        wigner = np.zeros((self.n,length*2+1))
        q = bit_bin*np.array([i for i in range(-length,length+1)])
        for num in range(self.n):
            for s in range(-photon_number,photon_number+1):
                for j in range(-length,length+1):
                    wigner[num,j] += np.exp(-(self.delta**2)/2*((2*s)**2)*np.pi)*np.exp(-(1/(2*self.delta**2)*((q[j]-2*s*np.sqrt(np.pi))**2)))
        plt.plot(wigner[0])
        plt.show()

    def __str__(self):
        return ""

    @property
    def n(self):
        return self._n