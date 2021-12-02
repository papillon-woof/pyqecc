from .stabilizer import *

class CONB(SC):
    def __init__(self,Codeinstances,Interleaver={},IID=True):
        self._IID = IID
        self.conc_length = len(Codeinstances) #連結符号の長さ
        self._Codeinstances = Codeinstances
        self._n = 0 #合計
        self._k = 0 #合計
        for c in self.Codeinstances:
            self._n += c.n
            self._k += c.k
        # Stabilizer
        self._H = np.zeros(self.n-self.k,2*self.n)
        nind = 0
        nkind = 0
        for c in self.Codeinstances:
            for h in c.H:
                self._H[nkind][nind:nind+c.n] = h[:nind]#X
                self._H[nkind][self.n+nind:self.n+nind+c.n] = h[nind:]#Z
            nind += c.n

        #2 ** kだけ必要なので，オーバーライドして減らす．
        def get_T(self,ind):
            '''
            indを受け取った時，以下をループし，足す．
            ind&(2 ** n - 1)
            ind = ind >> n
            '''
            T = np.zeros(2*self.n)
            nind = self.n
            #後ろから代入
            for i in reversed(range(self.code_length)):
                T_child = self.Codeinstances[i].get_T(ind&(2 ** self.Codeinstances[i].n - 1))
                T[nind-self.Codeinstances[i].n:nind] = T_child
                nind -= self.Codeinstances[i].n
                ind = ind>>self.Codeinstances[i].n #nだけシフトし，次の要素符号のTを取り出す準備
            return T

        def get_L(self,):
            L = np.zeros(2*self.n)
            nind = self.n
            #後ろから代入
            for i in reversed(range(self.code_length)):
                L_child = self.Codeinstances[i].get_L(ind&(2 ** self.Codeinstances[i].n - 1))
                L[nind-self.Codeinstances[i].n:nind] = L_child
                nind -= self.Codeinstances[i].n
                ind = ind>>self.Codeinstances[i].n #nだけシフトし，次の要素符号のTを取り出す準備
            return L

    @property
    def Codeinstances(self):
        return self._Codeinstances

    @property
    def conc_length(self):
        return self._conc_length

    @property
    def n_sum(self):
        return self._n_sum

    @property
    def k_sum(self):
        return self._k_sum

    @property
    def IID(self):
        return self._IID

#SのうちXだったらLxに拡張.e.g. XZIZX=LxLzILzLx\in 25ビット

class CONC(SC):
    def __init__(self,Codeinstances,Interleaver={},IID=True):
        '''
        array(SC): Codeinstances, dict: Interleaver, boolean: IID

        Codeinstances: 符号が格納．
        Codeinstancesは[[(1段目の符号インスタンス列)],[(2段目の符号インスタンス列)],[(3段目の符号インスタンス列)],...]
        と書きます．例えば，5量子ビットの連接符号の場合，
        例: Codeinstances = [[FIVE()],[FIVE() for i in range(5*1)],[FIVE() for i in range(5*1)]]
        と書きます．
        Interleaver: ビット間のインタリーバです．任意の量子ビットの配列で設定したいとき，インタリーブしてください．形式は転置行列で書きます．
        例: Interleaver = {1:T1,2:T2}
        '''
        self._IID = IID
        self._conc_length = len(Codeinstances) #連接符号の長さ
        self._Codeinstances = Codeinstances
        self._n_sum = [0]*self.conc_length #合計
        self._k_sum = [0]*self.conc_length #合計
        for l in range(self._conc_length):
            for c in self.Codeinstances[l]:
                self._n_sum[l] += c.n
                self._k_sum[l] += c.k
        for l in range(1,self._conc_length):
            if self._k_sum[l] != self._n_sum[l-1]:
                ValueError("Error:codelength mismatched")
        print(self.n_sum[self.conc_length-1],self.k_sum[0])
        super().__init__(self.n_sum[self.conc_length-1],self.k_sum[0])

    def calc_grobal_H(self):
        H = np.zeros(self.n-self.k,self.n)
        count = 0
        for l in range(self._conc_length):
            for c in self.Codeinstances[l]:
                H = np.zeros(self.n)
                for h in c.H:
                    for hi in h[:c.n]:
                        if 1==hi:
                            l+i
                        H[count]
                for k in range(self.Codeinstances[l].k):
                    self.Codeinstances[l].L[k]
                H.append()
    
    def BP_decode(self,s):
        for l in reversed(range(self.conc_length)):
            if self.IID:
                P_L = np.zeros((self.n_sum,4)) # 独立同一
            else:
                P_L = np.zeros((4 ** self.n_sum)) # 非独立同一
                ValueError("Not supported")
            qbit_position_count = 0
            for c in self.Codeinstances[l]:
                if l == self.conc_length-1:
                    c.set_P(self.P[qbit_position_count:qbit_position_count+c.n])
                else:
                    c.set_P(P_L[qbit_position_count:qbit_position_count+c.n])
                qbit_position_count+=c.n
            qbit_position_count = 0
            for c in self.Codeinstances[l]:
                L,P_L[qbit_position_count:qbit_position_count+c.k]= c.decoder(s,mode="ML",return_logical_error_probability=True)
        return L^T
    
    @property
    def Codeinstances(self):
        return self._Codeinstances

    @property
    def conc_length(self):
        return self._conc_length

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
    def __init__(self,Codeinstances,Interleaver={},IID=True):
        '''
        array(SC): Codeinstances, dict: Interleaver, boolean: IID

        Codeinstances: 符号が格納．
        Codeinstancesは[[(1段目の符号インスタンス列)],[(2段目の符号インスタンス列)],[(3段目の符号インスタンス列)],...]
        と書きます．例えば，5量子ビットの連接符号の場合，
        例: Codeinstances = [[FIVE()],[FIVE() for i in range(5*1)],[FIVE() for i in range(5*1)]]
        と書きます．
        Interleaver: ビット間のインタリーバです．任意の量子ビットの配列で設定したいとき，インタリーブしてください．形式は転置行列で書きます．
        例: Interleaver = {1:T1,2:T2}
        '''
        self._IID = IID
        self._conc_length = len(Codeinstances) #連接符号の長さ
        self._Codeinstances = Codeinstances
        self._n_sum = [0]*self.conc_length #合計
        self._k_sum = [0]*self.conc_length #合計
        for l in range(self._conc_length):
            for c in self.Codeinstances[l]:
                self._n_sum[l] += c.n
                self._k_sum[l] += c.k
        for l in range(1,self._conc_length):
            if self._k_sum[l] != self._n_sum[l-1]:
                ValueError("Error:codelength mismatched")
        print(self.n_sum[self.conc_length-1],self.k_sum[0])
        super().__init__(self.n_sum[self.conc_length-1],self.k_sum[0])

    def calc_grobal_H(self):
        H = np.zeros(self.n-self.k,self.n)
        count = 0
        for l in range(self._conc_length):
            for c in self.Codeinstances[l]:
                H = np.zeros(self.n)
                for h in c.H:
                    for hi in h[:c.n]:
                        if 1==hi:
                            l+i
                        H[count]
                for k in range(self.Codeinstances[l].k):
                    self.Codeinstances[l].L[k]
                H.append()
    
    def BP_decode(self,s):
        for l in reversed(range(self.conc_length)):
            if self.IID:
                P_L = np.zeros((self.n_sum,4)) # 独立同一
            else:
                P_L = np.zeros((4 ** self.n_sum)) # 非独立同一
                ValueError("Not supported")
            qbit_position_count = 0
            for c in self.Codeinstances[l]:
                if l == self.conc_length-1:
                    c.set_P(self.P[qbit_position_count:qbit_position_count+c.n])
                else:
                    c.set_P(P_L[qbit_position_count:qbit_position_count+c.n])
                qbit_position_count+=c.n
            qbit_position_count = 0
            for c in self.Codeinstances[l]:
                L,P_L[qbit_position_count:qbit_position_count+c.k]= c.decoder(s,mode="ML",return_logical_error_probability=True)
        return L^T
    
    @property
    def Codeinstances(self):
        return self._Codeinstances

    @property
    def conc_length(self):
        return self._conc_length

    @property
    def n_sum(self):
        return self._n_sum

    @property
    def k_sum(self):
        return self._k_sum

    @property
    def IID(self):
        return self._IID