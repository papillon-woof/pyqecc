import numpy as np
from .abstruct import *
from ..util import *
class SC(CODE):
    def __init__(self,n,k,H='random',T=None,L=None,P=None,iid=True,mode='HD'):
        super().__init__(n,k)
        self._mode = mode
        self._name = "stabilizer"
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

    def set_P(self,P,iid=True):
        if iid and P is not None:
            self._P=np.array([P.tolist()]*self.n)
        else:
            self._P = P

    def get_syndrome(self,e):
        return symplex_binary_inner_product(self._H,e)

    def get_T(self,ind):
        return self.T[arr2int(ind)]

    def get_S(self,ind):
        S = np.zeros(2*self.n,dtype='i1')
        for i in range(self.n-self.k):
            S+=(((ind>>i)&1)*self.H[i])
        return np.mod(S,2)

    def get_L(self,ind):
        L = np.zeros(2*self.n,dtype='i1')
        for i in range(2*self.k):
            L+=(((ind>>i)&1)*self.L[i])
        return np.mod(L,2)

    def in_S(self,b):
        return sum(gaussjordan(np.append(self.H,b).reshape(self.n-self.k+1,2*self.n))[self.n-self.k])==0

    #def _in_S(self,b):
    #    for i in range(self.n-self.k):
    #        if symplex_binary_inner_product(self.H[i],b)==1:
    #            return False
    #    return True

    def hard_decode(self,syndrome):
        T = self.get_T(syndrome)
        return T

    def ML_decode(self,syndrome):
        #decoding_metric: メトリック．受信語と通信路情報から計算する．
        T = self.get_T(syndrome)
        if self.n>self.ML_decoding_qubit_limit:
            ValueError("Error: The qubit n ="+str(self.n)+" is limited because of a large decoding complexity. You can change the qubit limit.")

        #Lについてビット全探索
        P_L = np.zeros(4**self.k)
        for li in range(4 ** self.k):
            SUM_P = np.zeros(4)
            L = self.get_L(li)
            #Sについてビット全探索
            for si in range(2 ** (self.n-self.k)):
                S = self.get_S(si)
                E = L^T^S
                Ptmp = 1
                for ei in range(self.n):
                    ind=E[ei]+2*E[ei+self.n]
                    Ptmp*=self.P[ei][ind]
                P_L[li]+=Ptmp
        l_ind = np.argmax(P_L)
        L = np.zeros(2*self.n,dtype='i1')
        for lj in range(2*self.k):
            L+=(((l_ind>>(lj))&1)*self.L[lj])

        return L^T

    def decode(self,syndrome):
        if self._mode=="ML":
            EE = self.ML_decode(syndrome)
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

    def __str__(self):
        output = ""
        output+="codename        :"+str(self.name)+"\n"
        output+="n               : "+str(self.n)+"\n"
        output+="k               : "+str(self.k)+"\n"
        output+="R               : "+str(self.R)+"\n"
        return output
