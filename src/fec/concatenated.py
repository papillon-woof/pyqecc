from .stabilizer import *
class CombCode(SC):
    '''
    符号を並列連結した符号を作る．
    code_instances: 符号インスタンスのlist.右から順に，通信路に近くなる．
    Interleaver: 量子ビット順序のインタリーバ．バースト誤り防止や順序を入れ替えたい場合に使用．
    IID:誤り分布が，独立同分布か．
    [110000|000000]
    [011000|000000]
    [000110|000000]
    [000011|000000]
    '''
    NAME = "COMBOLUTION_CODE"
    def __init__(self,code_instances,Interleaver={},P=None,iid=True,mode='ML',BITWISE = True):
        self.conc_length = len(code_instances) #連結符号の長さ
        self._code_instances = code_instances
        self._n = sum([c.n for c in self.code_instances])
        self._k = sum([c.k for c in self.code_instances])
        # Stabilizer
        self._H = np.zeros((self.n-self.k,2*self.n),dtype='i1')
        nind = 0
        nkind = 0
        for c in self.code_instances:
            for h in c.H:
                self._H[nkind][nind:nind+c.n] = h[:c.n]#X
                self._H[nkind][self.n+nind:self.n+nind+c.n] = h[c.n:2*c.n]#Z
                nkind += 1
            nind += c.n
        super().__init__(self.n,self.k,self._H,P=P,BITWISE = BITWISE,iid = iid,mode=mode)

    def get_T(self,beta):
        '''
        indを受け取った時，以下をループし，足す．
        ind&(2 ** n - 1)
        ind = ind >> n
        '''
        if len(beta) != self.n - self.k:
            raise ValueError("Length of beta is not matched number of stabilizer basis. Please check the length of beta.")
        T = np.zeros(2*self.n,dtype='i1')
        nind = self.n #シフトして，部分回復演算子を求めるため，後ろから求める．
        ind = self.n - self.k
        for c in reversed(self.code_instances):
            T_child = c.get_T(beta[ind-(c.n-c.k):ind]) #後ろの要素符号のインスタンスを取得
            T[nind-c.n:nind] = T_child[:c.n]#Xを代入
            T[nind-c.n+self.n:nind+self.n] = T_child[c.n:2*c.n]#Zを代入．
            nind -= c.n
            ind -= c.n-c.k #nだけシフトし，次の要素符号のTを取り出す準備
        return T

    #2 ** kだけ必要なので，オーバーライドして減らす．
    def get_L(self,alpha):
        L = np.zeros(2*self.n,dtype='i1')
        kind = 0
        nind = 0
        if type(alpha)==list:
            alpha = np.array(alpha)
        if type(alpha)==int or type(alpha)==np.int64:
            alpha = int2arr(alpha,2*self.k)
        if len(alpha) != 2 * self.k:
            raise ValueError("Length of alpha is not matched number of stabilizer basis. Please check the length of alpha.")
        for c in self.code_instances:
            L_child = c.get_L(np.concatenate([alpha[kind:kind+c.k],alpha[self.k + kind:self.k + kind+c.k]]))
            L[nind:nind+c.n] = L_child[:c.n]#Xを代入
            L[self.n + nind:self.n + nind+c.n] = L_child[c.n:2*c.n]#Zを代入
            kind += c.k
            nind += c.n
        return L

    def ML_decode(self,syndrome,**param):
        if "return_logical_error_probability" in param:
            self.return_logical_error_probability = param["return_logical_error_probability"]
        else:
            self.return_logical_error_probability = False
        L = np.zeros(X_AND_Z*self.n)
        T = self.get_T(syndrome)
        logical_error_probability = np.zeros(X_AND_Z*self.n,4)
        nkind = 0
        nind = 0
        for c in self.code_instances:
            l,p = c.ML_decode(syndrome[nkind:c.nk+nkind],**param)
            L[nind:nind+c.n],logical_error_probability[nind:nind+c.n][:] = l,blockwise_to_bitwise_probability(p)
            nkind += (c.n - c.k)
            nind += c.n
        if self.return_logical_error_probability:
            return L^T,blockwise_to_bitwise_probability(logical_error_probability)
        return L^T

    @property
    def code_instances(self):
        return self._code_instances

#SのうちXだったらLxに拡張.e.g. XZIZX=LxLzILzLx\in 25ビット

class ConcCode(SC):
    NAME = "CONCATENATED_CODE"
    def __init__(self,code_instances,interleaver={},P=None,iid=True,BITWISE = True):
        '''
        array(SC): code_instances, dict: Interleaver, boolean: IID
        code_instances: 符号が格納．
        code_instancesは[(1段目の符号インスタンス),(2段目の符号インスタンス列),...]
        と書く．例えば，符号２が[25,5]符号，符号１が[5,1]符号の場合，
        例: code_instances = [CODE1,CODE2]
        と書く．
        Interleaver: ビット間のインタリーバ．任意の量子ビットの配列で設定したいとき，インタリーブする．形式は転置行列で書く．
        例: Interleaver = {1:T1,2:T2}
        '''
        self._code_depth = len(code_instances) #連接符号の長さ
        self._code_instances = code_instances
        for d in range(0,self._code_depth-1):
            if self._code_instances[d+1].k != self._code_instances[d].n:
                raise ValueError("Error:codelength mismatched")
        self.H_depth = {}
        super().__init__(self._code_instances[self.conc_length-1].n,self._code_instances[0].k,H=self.calc_grobal_H(),P=P,iid = iid,BITWISE = BITWISE,mode = 'HD')

    def get_L(self,alpha):
        alpha = any2arr(alpha,2*self.k)
        if len(alpha) != 2 * self.k:
            raise ValueError("Length of alpha is not matched number of stabilizer basis. Please check the length of alpha.")
        L0 = self.code_instances[0].get_L(alpha)
        for d in range(self._code_depth-1):
            c0 = self.code_instances[d]
            c1 = self.code_instances[d+1]
            L1 = np.zeros(2*c1.n,dtype='i1')
            nind = 0
            for i in range(c0.n):
                for x_or_z in range(2):
                    if 1==L0[c0.n*x_or_z+i]:
                        Lind = np.zeros(2*c1.k,dtype='i1')
                        Lind[c1.k*x_or_z+i]=1
                        L1 += c1.get_L(Lind) #Lを格納．例えば，繰り返し符号ならL[0]=XXX，
            L0 = L1
        return L0

    def get_T(self,beta):
        beta0 = np.zeros(self.code_instances[0].n - self.code_instances[0].k,dtype='i1')
        for mother_ind,follower_ind in self.H_depth[0].items():
            beta0[follower_ind]=beta[mother_ind]
        T0 = self.code_instances[0].get_T(beta0)
        for d in range(self._code_depth-1):
            c0 = self.code_instances[d]
            c1 = self.code_instances[d+1]
            T1 = np.zeros(2*c1.n,dtype='i1')
            for i in range(c0.n):
                for x_or_z in range(2):
                    if 1==T0[c0.n*x_or_z+i]:
                        Lind = np.zeros(2*c1.k,dtype='i1')
                        Lind[c1.k*x_or_z+i]=1
                        T1 += c1.get_L(Lind) #Lを格納．例えば，繰り返し符号ならL[0]=XXX，
            beta1 = np.zeros(c1.n - c1.k,dtype='i1')
            for mother_ind,follower_ind in self.H_depth[d+1].items():
                beta1[follower_ind]=beta[mother_ind]
            T1 ^= c1.get_T(beta1)
            T0 = T1
        return T0

    def calc_grobal_H(self):
        H0 = self.code_instances[0].H
        for d in range(self._code_depth-1):
            c0 = self.code_instances[d]
            c1 = self.code_instances[d+1]
            H1 = np.zeros((c1.n-c1.k + H0.shape[0],2*c1.n),dtype='i1')
            nind = 0
            #for h0 in H0:
            for i in range(c0.n - c0.k):
                for bit in range(c0.n):
                    for x_or_z in range(2):
                        if 1==H0[i][c0.n*x_or_z+bit]:
                            Lind = np.zeros(2*c1.k,dtype='i1')
                            Lind[c1.k*x_or_z+bit]=1
                            H1[nind] += c1.get_L(Lind) #Lを格納．例えば，繰り返し符号ならL[0]=XXX，
                if d not in self.H_depth:
                    self.H_depth[d] = {}
                self.H_depth[d][nind] = i
                nind += 1 # 全体のHの行数をカウント
            for i in range(c1.n - c1.k):
                if d+1 not in self.H_depth:
                    self.H_depth[d+1] = {}
                self.H_depth[d+1][nind] = i
                H1[nind] = c1.H[i]
                nind+=1
            H0 = H1
        return H0

    def decode(self,syndrome,**param):
        self._mode = "BP"
        if "BP"==self.mode:
            EE=self.BP_decode(syndrome)
        else:
            EE = super().decode(syndrome,**param)
        return EE

    def BP_decode(self,syndrome):
        '''
        BP_decoding:
        input: syndrome
        output: L^T
        '''

        # Estimation for Logical error
        P = self.bitwise_p
        for d in reversed(range(self._code_depth)):
            c1 = self.code_instances[d]
            beta1 = np.zeros(c1.n - c1.k,dtype='i1')
            for mother_ind,follower_ind in self.H_depth[d].items():
                beta1[follower_ind]=syndrome[mother_ind]
            c1.set_error_probability(P)
            L0,P = c1.decode(beta1,mode="ML",return_logical_error_probability=True)
        for d in range(self._code_depth-1):
            c0 = self.code_instances[d]
            c1 = self.code_instances[d+1]
            L1 = np.zeros(2*c1.n,dtype='i1')
            for i in range(c0.n):
                for x_or_z in range(2):
                    if 1==L0[c0.n*x_or_z+i]:
                        Lind = np.zeros(2*c1.k,dtype='i1')
                        Lind[c1.k*x_or_z+i]=1
                        L1 += c1.get_L(Lind) #Lを格納．例えば，繰り返し符号ならL[0]=XXX，
            L0 = L1
        return L0^self.get_T(syndrome)

    @property
    def code_instances(self):
        return self._code_instances

    @property
    def conc_length(self):
        return self._code_depth

    @property
    def iid(self):
        return self._iid