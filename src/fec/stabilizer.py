import numpy as np
from .abstruct import *
from ..util import *

class SC(CODE):
    '''
    '''
    NAME = "stabilizer code"
    USAGE = "Text"
    ML_DECODING_QUBITS_LIMIT = 15
    def __init__(self,n,k,H='random',T=None,L=None,P=None,iid=True,mode='HD'):
        super().__init__(n,k)
        self._mode = mode
        self._H = H.astype('i1')
        self._T = T
        self._L = L
        self._P = self.set_P(P)
        self._LUT = {}
        #self.enc_circuit = None
        #self.dec_circuit = None

    def set_P(self,P,iid=True):
        if iid and P is not None:
            self._P=np.array([P.tolist()]*self.n)
        else:
            self._P = P

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
        if type(alpha)==list:
            alpha = np.array(alpha)
        if type(alpha)==int or type(alpha)==np.int64:
            alpha = int2arr(alpha,2*self.k)
        if len(alpha) != 2 * self.k:
            raise ValueError("Length of alpha is not matched number of stabilizer basis. Please check the length of alpha.")
        L = np.zeros(2*self.n,dtype='i1')
        #LX or LZ
        for i in range(2):
            for j in range(self.k):
                L+=alpha[i*self.k+j]*self.L[j][i]
        return np.mod(L,2)

    def in_S(self,b):
        return sum(gaussjordan(np.append(self.H,b).reshape(self.n-self.k+1,2*self.n))[self.n-self.k])==0

    def set_LUT(self):
        for i in range(2 ** (self.n - self.k)):
            self._LUT[i] = self.ML_decode(int2arr(i,(self.n - self.k)))

    def hard_decode(self,syndrome):
        T = self.get_T(syndrome)
        return T.astype("i1")

    def ML_decode(self,syndrome,**param):
        #decoding_metric: メトリック．受信語と通信路情報から計算する．
        if "return_logical_error_probability" in param:
            self.return_logical_error_probability = param["return_logical_error_probability"]
        else:
            self.return_logical_error_probability = False
        T = self.get_T(syndrome)
        if self.n>self.ML_DECODING_QUBITS_LIMIT:
            ValueError("Error: The qubit n ="+str(self.n)+" is limited because of a large decoding complexity. You can change the qubit limit.")
        #Lについてビット全探索
        P_L = np.zeros(2**(2*self.k))
        for li in range(2 ** (2*self.k)):
            L = self.get_L(li)
            #Sについてビット全探索
            for si in range(2 ** (self.n-self.k)):
                S = self.get_S(si)
                E = L^T^S
                Ptmp = 1
                for ei in range(self.n):
                    ind=E[ei]+(E[ei+self.n]<<1)
                    Ptmp*=self.P[ei][ind]
                P_L[li]+=Ptmp
        l_ind = np.argmax(P_L)
        L = self.get_L(l_ind)
        if self.return_logical_error_probability:
            return L^T,P_L
        return L^T

    def LUT_decode(self,syndrome):
        return self.LUT[arr2int(syndrome)]

    def decode(self,syndrome,**param):
        if self._mode=="ML":
            EE = self.ML_decode(syndrome)
        if self._mode=="ML_LUT":
            if self.LUT == {}:
                self.set_LUT()
            EE = self.LUT_decode(syndrome)
        elif self._mode=="HD":
            EE = self.hard_decode(syndrome)
        return EE

    @property
    def P(self):
        return self._P

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

    def __str__(self):
        output = ""
        output+="NAME            :"+str(self.name)+"\n"
        output+="N               : "+str(self.n)+"\n"
        output+="K               : "+str(self.k)+"\n"
        output+="R               : "+str(self.R)+"\n"
        return output
