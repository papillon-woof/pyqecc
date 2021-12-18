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
    def __init__(self,code_instances,Interleaver={},iid=True):
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
        self._H = np.zeros((self.n-self.k,2*self.n))
        nind = 0
        nkind = 0
        for c in self.code_instances:
            for h in c.H:
                self._H[nkind][nind:nind+c.n] = h[:c.n]#X
                self._H[nkind][self.n+nind:self.n+nind+c.n] = h[c.n:2*c.n]#Z
                nkind += 1
            nind += c.n
        print(self._H)    

    def get_T(self,beta):
        '''
        indを受け取った時，以下をループし，足す．
        ind&(2 ** n - 1)
        ind = ind >> n
        '''
        if len(beta) != self.n - self.k:
            raise ValueError("Length of beta is not matched number of stabilizer basis. Please check the length of beta.")
        T = np.zeros(2*self.n)
        nind = self.n #シフトして，部分回復演算子を求めるため，後ろから求める．
        ind = 0
        for c in reversed(self.code_instances):
            print(self.n,nind,nind-2*c.n,nind-c.n)
            T_child = c.get_T(beta[ind:ind+c.n-c.k]) #後ろの要素符号のインスタンスを取得
            T[nind-c.n:nind] = T_child[:c.n]#Xを代入
            T[nind-c.n+self.n:nind+self.n] = T_child[c.n:2*c.n]#Zを代入．
            nind -= c.n
            ind = c.n-c.k #nだけシフトし，次の要素符号のTを取り出す準備
        return T
    
    # #2 ** kだけ必要なので，オーバーライドして減らす．
    # def get_L(self,ind):
    #     L = np.zeros(2*self.n)
    #     nind = self.n
    #     #後ろから代入
    #     for c in reversed(self.code_instances):
    #         L_child = c.get_L(ind&(2 ** c.n - 1))
    #         L[nind-2*c.n:nind-c.n] = L_child[:c.n]#Xを代入
    #         L[nind-c.n:nind] = L_child[c.n:2*c.n]#Zを代入
    #         nind -= c.n
    #         ind = ind>>c.n #nだけシフトし，次の要素符号のTを取り出す準備
    #     return L

    #2 ** kだけ必要なので，オーバーライドして減らす．
    def get_L(self,alpha):
        L = np.zeros(2*self.n)
        nkind = self.n - self.k
        nind = self.n
        if len(alpha) != 2 * (self.n - self.k):
            raise ValueError("Length of alpha is not matched number of stabilizer basis. Please check the length of alpha.")
        #後ろから代入
        for c in self.code_instances:
            L_child = c.get_L(np.c_[alpha[nkind:nkind+c.n-c.k],alpha[self.n - self.k + nkind:self.n - self.k + nkind+c.n-c.k]])
            L[nind:nind+c.n] = L_child[:c.n]#Xを代入
            L[self.n + nkind:self.n + nkind+c.n] = L_child[c.n:2*c.n]#Zを代入
            nkind += c.n - c.k
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
        for l in range(1,self._code_depth):
            if self._k_sum[l] != self._n_sum[l-1]:
                ValueError("Error:codelength mismatched")
        print(self.n_sum[self.conc_length-1],self.k_sum[0])
        super().__init__(self.n_sum[self.conc_length-1],self.k_sum[0])

    def calc_grobal_H(self):
        H = np.zeros(self.n-self.k,self.n)
        count = 0
        for c in self.code_instances:
            for h in c.H:
                for hi in h[:c.n]:
                    if 1==hi:
                        H[count] = c.L[hi] #Lを格納．例えば，繰り返し符号ならL[0]=XXX，
                    H[count]
                count += 1 # 全体のHの行数をカウント
            for k in range(self.code_instances[l].k):
                self.code_instances[l].L[k]
            H.append()

    def BP_decode(self,s):
        for l in reversed(range(self.conc_length)):
            if self.IID:
                P_L = np.zeros((self.n_sum,4)) # 独立同一
            else:
                P_L = np.zeros((4 ** self.n_sum)) # 非独立同一
                ValueError("Not supported")
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
    def n_sum(self):
        return self._n_sum

    @property
    def k_sum(self):
        return self._k_sum

    @property
    def IID(self):
        return self._IID

class CONC(SC):
    NAME = "Concatenated Code"
    USAGE = ""
    def __init__(self,code_instances,Interleaver={},IID=True):
        '''
        array(SC): code_instances, dict: Interleaver, boolean: IID
        code_instances: 符号が格納．
        code_instancesは[[(1段目の符号インスタンス列)],[(2段目の符号インスタンス列)],[(3段目の符号インスタンス列)],...]
        と書きます．例えば，5量子ビットの連接符号の場合，
        例: code_instances = [[FIVE()],[FIVE() for i in range(5*1)],[FIVE() for i in range(5*1)]]
        と書きます．
        Interleaver: ビット間のインタリーバです．任意の量子ビットの配列で設定したいとき，インタリーブしてください．形式は転置行列で書きます．
        例: Interleaver = {1:T1,2:T2}
        '''
        self._iid = iid
        self._code_depth = len(code_instances) #連接符号の長さ
        self._code_instances = code_instances
        self._n_sum = [0]*self.conc_length #合計
        self._k_sum = [0]*self.conc_length #合計
        for l in range(self._code_depth):
            for c in self.code_instances[l]:
                self._n_sum[l] += c.n
                self._k_sum[l] += c.k
        for l in range(1,self._code_depth):
            if self._k_sum[l] != self._n_sum[l-1]:
                ValueError("Error:codelength mismatched")
        print(self.n_sum[self.conc_length-1],self.k_sum[0])
        super().__init__(self.n_sum[self.conc_length-1],self.k_sum[0])

    def calc_grobal_H(self):
        H = np.zeros(self.n-self.k,self.n)
        count = 0
        for l in range(self._code_depth):
            for c in self.code_instances[l]:
                H = np.zeros(self.n)
                for h in c.H:
                    for hi in h[:c.n]:
                        if 1==hi:
                            l+i
                        H[count]
                for k in range(self.code_instances[l].k):
                    self.code_instances[l].L[k]
                H.append()
    
    def BP_decode(self,s):
        for l in reversed(range(self.conc_length)):
            if self.IID:
                P_L = np.zeros((self.n_sum,4)) # 独立同一
            else:
                P_L = np.zeros((4 ** self.n_sum)) # 非独立同一
                ValueError("Not supported")
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
    def n_sum(self):
        return self._n_sum

    @property
    def k_sum(self):
        return self._k_sum

    @property
    def IID(self):
        return self._IID