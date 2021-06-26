import numpy as np
from .abstruct import *
from ..util import *
class SC(CODE):
    def __init__(self,n,k,H='random',T=None,L=None,P=None,iid=True):
        self._name = "stabilizer"
        self._n = n
        self._k = k
        self._R = self._k/self._n
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


    #量子情報ビット
    def get_enc_circ(self):
        return self.enc_circuit

    #シンドロームと確率分布を受け取る．
    def get_dec_circ(self):
        return self.dec_circuit

    # hard decision
    def get_syndrome(self,e):
        return symplex_binary_inner_product(self._H,e)

    def get_T(self,ind):
        T = np.zeros(2*self.n,dtype='i1')
        for i in range(self.n - self.k):
            T+=(ind[i]*self.T[i])
        return np.mod(T,2)

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
        #print(np.append(self.H,b).reshape(self.n-self.k+1,2*self.n))
        #print(gaussjordan(np.append(self.H,b).reshape(self.n-self.k+1,2*self.n)))
        return sum(gaussjordan(np.append(self.H,b).reshape(self.n-self.k+1,2*self.n))[self.n-self.k])==0

    def hard_decode(self,syndrome):
        T = self.get_T(syndrome)
        return T

    def ML_decode(self,syndrome):
        #decoding_metric: メトリック．受信語と通信路情報から計算する．
        T = self.get_T(syndrome)
        #print(111111,syndrome,T,self.get_syndrome(T))
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
                E = T^S^L
                #if 0==sum(E - np.array([0,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0,0])):
                #    print(T,S,L)
                #E = ()のうち，確率を入力
                Ptmp = 1
                for ei in range(self.n):
                    ind=E[ei]+2*E[ei+self.n]
                    Ptmp*=self.P[ei][ind]
                P_L[li]+=Ptmp
        l_ind = np.argmax(P_L)
        ##print(P_L)
        L = np.zeros(2*self.n,dtype='i1')
        for lj in range(2*self.k):
            L+=(((l_ind>>(lj))&1)*self.L[lj])

        #print("L,LT,T",self.get_syndrome(L),self.get_syndrome(L^T),self.get_syndrome(T))
        return L^T

    def decode(self,syndrome,mode='HD'):
        if mode=="ML":
            EE = self.ML_decode(syndrome)
        elif mode=="HD":
            EE = self.hard_decode(syndrome)
        return EE

    @property
    def P(self):
        return self._P

    def __str__(self):
        output = ""
        output+="codename        :"+str(self.name)+"\n"
        output+="n               : "+str(self.n)+"\n"
        output+="k               : "+str(self.k)+"\n"
        output+="R               : "+str(self.R)+"\n"
        return output
