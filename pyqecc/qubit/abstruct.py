from abc import ABCMeta, abstractmethod
import numpy as np

# measure
Measure = None
SW_gate = None
# 1 qubit gate
M_gate = np.array([[1,0],[0,0]],dtype=np.complex128)
I_gate = np.array([[1,0],[0,1]],dtype=np.complex128)
X_gate = np.array([[0,1],[1,0]],dtype=np.complex128)
Y_gate = np.array([[0,1j],[-1j,0]],dtype=np.complex128)
H_gate = np.array([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]],dtype=np.complex128)
Z_gate = np.array([[1,0],[0,-1]],dtype=np.complex128)
T_gate = np.array([[1,0],[0,np.exp(1j*np.pi/4)]],dtype=np.complex128)
# 2 qubit gate
CX_gate = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],dtype=np.complex128)
CZ_gate = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]],dtype=np.complex128)

# bosonic qubit
class Qubits(metaclass=ABCMeta):
    _name = "qubits"
    gate_fundamental_dictionary = {
        "M":M_gate,
        "X":X_gate,
        "Y":Y_gate,
        "Z":Z_gate,
        "T":T_gate,
        "CX":CX_gate,
        "CZ":CZ_gate,
        "SW":SW_gate,
        "I":I_gate,
    }
    def __init__(self,n):
        self.n = n
        self.qubits = np.zeros(2 ** n,dtype=np.complex128);self.qubits[0] = 1 # |000...000>準備
        self.gate_dictionary = {}
        # 動作用の行列を準備
        for k,gate in self.gate_fundamental_dictionary.items():
            if "SW"==k:
                continue
            gate_tmp = gate.copy()
            for idx in range(self.n-1):
                gate_tmp = np.kron(gate_tmp,self.gate_fundamental_dictionary["I"])
            self.gate_dictionary[k] = gate_tmp

    def M(self,first_idx=-1,second_idx=-1):
        r = np.random.rand()
        self.SW(first_idx=first_idx,second_idx=second_idx)
        p = np.real(np.dot(self.qubits.conjugate(),np.dot(self.gate_dictionary["M"],self.qubits)))
        self.SW(first_idx=first_idx,second_idx=second_idx)
        if r<p:
            return 0
        else:
            return 1

    def SW(self,first_idx=-1, second_idx=-1, therd_idx=-1, fourth_idx=-1):
        rep = {}
        if first_idx == second_idx:
            return 1.0
        for i in range(2 ** self.n):
            repBits = i
            repBits = repBits^(((repBits>>(self.n-1-second_idx))&1)<<(self.n-1-first_idx))
            repBits = repBits^(((repBits>>(self.n-1-first_idx))&1)<<(self.n-1-second_idx))
            repBits = repBits^(((repBits>>(self.n-1-second_idx))&1)<<(self.n-1-first_idx))
            if i>repBits:
                minTemp = i
                maxTemp = repBits
            else:
                maxTemp = i
                minTemp = repBits
            if minTemp not in rep and minTemp != maxTemp:
                rep[minTemp]=maxTemp
        for k,v in rep.items():
            self.qubits[k],self.qubits[v] = self.qubits[v],self.qubits[k]
        return 

    def gate(self,name,first_idx=-1,second_idx=-1):
        if "SW" == name:
            return self.SW(first_idx=first_idx,second_idx=second_idx)
        elif "M" == name:
            return self.M(first_idx=0,second_idx=first_idx)
        elif len(name) == 1:
            if 0 > first_idx or first_idx >= self.n:
                raise ValueError("invalid first index")
            if name not in self.gate_dictionary.keys():
                raise ValueError("invalid name")
            self.SW(first_idx=first_idx,second_idx=0)
            self.qubits = np.dot(self.gate_dictionary[name],self.qubits)
            self.SW(first_idx=first_idx,second_idx=0)
        
        elif len(name) == 2:
            if (0 > second_idx or second_idx >= self.n) and len(name) >= 2:
                raise ValueError("invalid second index")
            self.SW(first_idx=first_idx,second_idx=0,therd_idx=second_idx, fourth_idx=1)
            self.gate_dictionary[name]
            self.SW(first_idx=first_idx,second_idx=0,therd_idx=second_idx, fourth_idx=1)
        return 0

    def __str__(self,ket=True):
        if ket:
            output = ""
            for i in range(2 ** self.n):
                if self.qubits[i]!=0:
                    output+= (("(" if self.qubits[i].imag!=0 and self.qubits[i].real!=0 else "")+str(self.qubits[i].real)+"+" if self.qubits[i].imag>=0 else "") + str(self.qubits[i].imag)+"j"+("(" if self.qubits[i].imag!=0 and self.qubits[i].real!=0 else "")+"|"+str(bin(i))[2:].zfill(self.n)+">"+"+"
            output = output[:-1]
        return output
    