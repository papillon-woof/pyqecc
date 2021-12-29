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
    def __init__(self,code_instances,Interleaver={},iid=True,mode='ML'):
        self._mode = mode
        self._L = None
        self._iid = iid
        self.conc_length = len(code_instances) #連結符号の長さ
        self._code_instances = code_instances
        self._n = 0 #合計
        self._k = 0 #合計
        for c in self.code_instances:
            self._n += c.n
            self._k += c.k
        self._R = self.k/self.n
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
        #後ろから代入
        for c in self.code_instances:
            L_child = c.get_L(np.concatenate([alpha[kind:kind+c.k],alpha[self.k + kind:self.k + kind+c.k]]))
            L[nind:nind+c.n] = L_child[:c.n]#Xを代入
            L[self.n + nind:self.n + nind+c.n] = L_child[c.n:2*c.n]#Zを代入
            kind += c.k
            nind += c.n
        return L

    @property
    def code_instances(self):
        return self._code_instances

    @property
    def code_num(self):
        return self._code_num

    @property
    def n_sum(self):
        return self._n

    @property
    def k_sum(self):
        return self._k

    @property
    def iid(self):
        return self._iid
    

#SのうちXだったらLxに拡張.e.g. XZIZX=LxLzILzLx\in 25ビット

class ConcCode(SC):
    NAME = "CONCATENATED_CODE"
    def __init__(self,code_instances,interleaver={},iid=True):
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
        self._iid = iid # i.i.dかどうか．
        self._code_depth = len(code_instances) #連接符号の長さ
        self._code_instances = code_instances
        for d in range(0,self._code_depth-1):
            if self._code_instances[d+1].k != self._code_instances[d].n:
                ValueError("Error:codelength mismatched")
        super().__init__(self._code_instances[self.conc_length-1].n,self._code_instances[0].k)
        self.H_depth = {}
        self._H = self.calc_grobal_H()

    def get_L(self,alpha):
        if type(alpha)==list:
            alpha = np.array(alpha)
        if type(alpha)==int or type(alpha)==np.int64:
            alpha = int2arr(alpha,2*self.k)
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
        super().decode(syndrome,syndrome,**param)
        if "BP"==self.mode:
            EE=self.BP_decode(syndrome)
        return EE

    def BP_decode(self,s):
        for l in reversed(range(self.conc_length)):
            if self.IID:
                P_L = np.zeros((self.n_sum,4)) # 独立同一
            else:
                ValueError("N/A: non-i.i.d")
                P_L = np.zeros((4 ** self.n_sum)) # 非独立同一
            qbit_position_count = 0
            for c in self.code_instances[l]:
                if l == self.conc_length-1:
                    c.set_P(self.P[qbit_position_count:qbit_position_count+c.n])
                else:
                    c.set_P(P_L[qbit_position_count:qbit_position_count+c.n])
                qbit_position_count+=c.n
            qbit_position_count = 0
            for c in self.code_instances[l]:
                L,P_L[qbit_position_count:qbit_position_count+c.k]= c.decoder(s,mode="ML",return_logical_error_probability=True)
        return L^T

    @property
    def code_instances(self):
        return self._code_instances

    @property
    def conc_length(self):
        return self._code_depth

    @property
    def IID(self):
        return self._IID