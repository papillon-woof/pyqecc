import numpy as np
from .abstruct import *
from ..util import *

class SC(CODE):
    '''
    '''
    NAME = "stabilizer code"
    USAGE = "Text"
    ML_DECODING_QUBITS_LIMIT = 10
    decoder_output = {
            "LT" : None,
            "L" : None,
            "T" : None,
            "LOGICAL_ERROR_PROBABILITY" : None,
    }
    def __init__(
        self,
        n,
        k,
        H = None,
        T = None,
        L = None,
        P = None,
        mode = 'HD',
        BITWISE = True,
        ):
        super().__init__(n,k)
        self._mode = mode
        self._H = H
        self._T = T
        self._L = L
        self.BITWISE = BITWISE
        self.set_error_probability(P,self.BITWISE)
        self._LUT = {}
        #self.enc_circuit = None
        #self.dec_circuit = None

    def set_error_probability(self,P,BITWISE=True,iid=True):
        if not P is None:
            if BITWISE:
                if iid:
                    P = np.array([P for i in range(self.n)])
                self.set_bitwise_p(P)
                if self.ML_DECODING_QUBITS_LIMIT<self.n:
                    print("The num. of qubit exceed ML_DECODING_QUBITS_LIMIT.")
                    self._blockwise_p = False
                else:
                    self.set_blockwise_p(bitwise_to_blockwise_probability(P)) # By approximate the bitwise probability.
            else:
                self.set_blockwise_p(P)
                self.set_bitwise_p(blockwise_to_bitwise_probability(P))
        else:
            print("The input error probability is eligable. Please chack the input error probability.")

    def set_blockwise_p(self,blockwise_p):
        self._blockwise_p = blockwise_p

    def set_bitwise_p(self,bitwise_p):
        self._bitwise_p = bitwise_p

    def get_syndrome(self,e):
        return symplex_binary_inner_product(self._H,e)

    def get_T(self,ind):
        return self.T[arr2int(ind)] #LUTでの計算? BPでの計算もあり?LDPCについて学ぶ．

    def get_S(self,ind_list):
        S = np.zeros(2*self.n,dtype='i1')
        if type(ind_list) == int:
            ind_list = int2arr(ind_list,self.n - self.k)
        for i in range(self.n-self.k):
            S+=ind_list[i]*self.H[i]
        return np.mod(S,2)

    def get_L(self,alpha):
        '''
        alpha
        [LX1,LX2,LX3,...,LX(n-k)|LZ1,LZ2,LZ3,...,LZ(n-k)]
        '''
        alpha = any2arr(alpha,X_OR_Z*self.k)
        L = np.zeros(2*self.n,dtype='i1')
        #LX or LZ
        for j in range(self.k):
            for i in range(X_OR_Z):
                L^=alpha[i*self.k+j]*self.L[j][i]
        return np.mod(L,2)

    def in_S(self,b):
        return sum(gaussjordan(np.append(self.H,b).reshape(self.n-self.k+1,2*self.n))[self.n-self.k])==0

    def set_LUT(self):
        for i in range(2 ** (self.n - self.k)):
            self._LUT[i] = self.ML_decode(int2arr(i,(self.n - self.k)))

    def get_error_probability(self,E):
        return self.blockwise_p[arr2int(E)]

    def ML_decode(self,syndrome,**param): #logical error の複数ビット対応が必要．
        if self.n>self.ML_DECODING_QUBITS_LIMIT:
            raise ValueError("Error: The qubit n ="+str(self.n)+" is limited because of a large decoding complexity. You can change the qubit limit.")
        T = self.get_T(syndrome)
        logical_error_probability = np.zeros(2**(2*self.k))
        #L,Sについてビット全探索
        for lind in range(2 ** (2*self.k)):
            L = self.get_L(lind)
            for sind in range(2 ** (self.n-self.k)):
                S = self.get_S(sind)
                E = L^T^S
                logical_error_probability[lind]+=self.get_error_probability(E)
        llind = np.argmax(logical_error_probability)
        self.decoder_output["L"] = self.get_L(llind)
        self.decoder_output["T"] = self.get_T(syndrome)
        self.decoder_output["LT"] = self.decoder_output["L"]^self.decoder_output["T"]
        self.decoder_output["LOGICAL_ERROR_PROBABILITY"] = logical_error_probability

    def hard_decode(self,syndrome):
        self.decoder_output["LT"] = self.get_T(syndrome)

    def LUT_decode(self,syndrome):
        self.decoder_output["LT"] = self.LUT[arr2int(syndrome)]

    def decode(self,syndrome,mode=False,**param,):
        if not mode is False:
            self._mode = mode
        if self._mode=="ML":
            self.ML_decode(syndrome,**param)
        if self._mode=="ML_LUT":
            if self.LUT == {}:
                self.set_LUT()
            self.LUT_decode(syndrome,**param)
        elif self._mode=="HD":
            self.hard_decode(syndrome,**param)
        return self.decoder_output

    @property
    def blockwise_p(self):
        return self._blockwise_p

    @property
    def bitwise_p(self):
        return self._bitwise_p

    @property
    def L(self):
        return self._L

    @property
    def T(self):
        return self._T

    @property
    def H(self):
        return self._H

    @property
    def LUT(self):
        return self._LUT

    @property
    def mode(self):
        return self._mode

    def __str__(self):
        output = ""
        output+="NAME            :"+str(self.name)+"\n"
        output+="N               : "+str(self.n)+"\n"
        output+="K               : "+str(self.k)+"\n"
        output+="R               : "+str(self.R)+"\n"
        output+="DECODING_MODE   : "+str(self.mode)+"\n"
        return output
