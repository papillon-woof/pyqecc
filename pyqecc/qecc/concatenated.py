from .stabilizer import *


class ParaCode(SC):
    """
    Generate a paralellized code
    code_instances: a list of code instance which is near the channel from right to left.
    e.g. bit-flip code with parallelized number of two.
    [110000|000000]
    [011000|000000]
    [000110|000000]
    [000011|000000]
    """

    NAME = "COMBOLUTION_CODE"

    def __init__(
        self, code_instances, P=None, mode="ML", BITWISE=True
    ):  # Interleaver={}, 追加予定．
        self._num_code = len(code_instances)  # 連結符号の長さ
        self._code_instances = code_instances  # 要素符号インスタンス
        super().__init__(
            sum([c.n for c in self.code_instances]),
            sum([c.k for c in self.code_instances]),
            H=None,
            P=P,
            BITWISE=BITWISE,
            mode=mode,
        )
        # Stabilizer
        nind, nkind = 0, 0
        self._H = np.zeros((self.nk, X_OR_Z * self.n), dtype="i1")  # パリティ検査行列
        for c in self.code_instances:
            for h in c.H:
                self._H[nkind][nind : nind + c.n] = h[: c.n]  # X
                self._H[nkind][self.n + nind : self.n + nind + c.n] = h[
                    c.n : X_OR_Z * c.n
                ]  # Z
                nkind += 1
            nind += c.n
        if not BITWISE:
            raise ValueError("No implimentation. Please wait to update.")

    def get_T(self, beta):
        if len(beta) != self.nk:
            raise ValueError(
                "Length of beta is not matched number of stabilizer basis. Please check the length of beta."
            )
        T = np.zeros(X_OR_Z * self.n, dtype="i1")
        nind = self.n  # シフトして，部分回復演算子を求めるため，後ろから求める．
        nkind = self.nk
        for c in reversed(self.code_instances):
            T_child = c.get_T(beta[nkind - (c.nk) : nkind])  # 後ろの要素符号のインスタンスを取得
            T[nind - c.n : nind] = T_child[: c.n]  # Xを代入
            T[nind - c.n + self.n : nind + self.n] = T_child[
                c.n : X_OR_Z * c.n
            ]  # Zを代入．
            nind, nkind = nind - c.n, nkind - c.nk
        return T

    def get_L(self, alpha):
        L = np.zeros(X_OR_Z * self.n, dtype="i1")
        kind, nind = 0, 0
        alpha = any2arr(alpha, X_OR_Z * self.k)
        for c in self.code_instances:
            L_child = c.get_L(
                np.concatenate(
                    [
                        alpha[kind : kind + c.k],
                        alpha[self.k + kind : self.k + kind + c.k],
                    ]
                )
            )
            L[nind : nind + c.n] = L_child[: c.n]  # Xを代入
            L[self.n + nind : self.n + nind + c.n] = L_child[c.n : X_OR_Z * c.n]  # Zを代入
            kind, nind = kind + c.k, nind + c.n
        return L

    def ML_decode(self, syndrome, **param):
        L = np.zeros(X_OR_Z * self.n, dtype="i1")
        logical_error_probability = np.zeros((self.k, 4))
        nkind, nind, kind = 0, 0, 0
        for c in self.code_instances:
            c.set_channel_param(
                self.bitwise_error_probability[nind : nind + c.n], iid=False
            )
            c_decoder_output = c.decode(syndrome[nkind : c.n - c.k + nkind], mode="ML")
            L[X_OR_Z * nind : X_OR_Z * nind + X_OR_Z * c.n] = c_decoder_output["L"]
            c_logical_error_probability = c_decoder_output["LOGICAL_ERROR_PROBABILITY"]
            logical_error_probability[kind : kind + c.k][
                :
            ] = blockwise_to_bitwise_error_probability(c_logical_error_probability)
            nind, kind, nkind = c.n + nind, c.k + kind, c.nk + nkind
        self.decoder_output["L"] = L
        self.decoder_output["T"] = self.get_T(syndrome)
        self.decoder_output["LT"] = self.decoder_output["L"] ^ self.decoder_output["T"]
        if self.BITWISE:
            self.decoder_output[
                "BIT_WISE_LOGICAL_ERROR_PROBABILITY"
            ] = logical_error_probability
        else:
            self.decoder_output[
                "LOGICAL_ERROR_PROBABILITY"
            ] = bitwise_to_blockwise_error_probability(logical_error_probability)

    @property
    def code_instances(self):
        return self._code_instances

    @property
    def num_code(self):
        return self._num_code


# SのうちXだったらLxに拡張.e.g. XZIZX=LxLzILzLx\in 25ビット
class ConcCode(SC):
    def __init__(
        self, code_instances, P=None, BITWISE=True, mode="BP"
    ):  # ,interleaver={}
        """
        array(SC): code_instances, dict: Interleaver, boolean: IID
        code_instances: 符号が格納．
        code_instancesは[(1段目の符号インスタンス),(2段目の符号インスタンス列),...]
        と書く．例えば，符号２が[25,5]符号，符号１が[5,1]符号の場合，
        例: code_instances = [CODE1,CODE2]
        と書く．
        Interleaver: ビット間のインタリーバ．任意の量子ビットの配列で設定したいとき，インタリーブする．形式は転置行列で書く．
        例: Interleaver = {1:T1,2:T2}
        """
        self.NAME = "CONCATENATED_CODE"
        self._code_depth = len(code_instances)  # 連接符号の長さ
        self._code_instances = code_instances
        for d in range(0, self._code_depth - 1):
            if self.code_instances[d + 1].k != self.code_instances[d].n:
                raise ValueError("Error:codelength mismatched")
        self.H_depth = {}
        super().__init__(
            self.code_instances[self.code_depth - 1].n,
            self.code_instances[0].k,
            H=self.calc_grobal_H(),
            P=P,
            BITWISE=BITWISE,
            mode=mode,
        )

    def get_mother_operator(self, U, depth_U):
        # depth_U : depth
        # U : operator
        if depth_U == self.code_depth:
            return U
        U0 = U
        for d in range(depth_U, self.code_depth - 1):
            c0 = self.code_instances[d]
            c1 = self.code_instances[d + 1]
            U1 = np.zeros(X_OR_Z * c1.n, dtype="i1")
            for i in range(c0.n):
                for x_or_z in range(X_OR_Z):
                    if 1 == U0[c0.n * x_or_z + i]:
                        Lind = np.zeros(X_OR_Z * c1.k, dtype="i1")
                        Lind[c1.k * x_or_z + i] = 1
                        U1 ^= c1.get_L(Lind)  # Lを格納．例えば，繰り返し符号ならL[0]=XXX，
            U0 = U1
        return U0

    def get_L(self, alpha):
        alpha = any2arr(alpha, X_OR_Z * self.k)
        if len(alpha) != X_OR_Z * self.k:
            raise ValueError(
                "Length of alpha is not matched number of stabilizer basis. Please check the length of alpha."
            )
        L0 = self.code_instances[0].get_L(alpha)
        return self.get_mother_operator(L0, 0)

    def get_T(self, beta):
        T = np.zeros(X_OR_Z * self.n, dtype="i1")
        for d in range(self._code_depth):
            beta0 = np.zeros(self.code_instances[d].nk, dtype="i1")
            for mother_ind, follower_ind in self.H_depth[d].items():
                beta0[follower_ind] = beta[mother_ind]
            T ^= self.get_mother_operator(self.code_instances[d].get_T(beta0), d)
        return T

    def calc_grobal_H(self):
        nind = 0
        n = self._code_instances[self.code_depth - 1].n
        k = self._code_instances[0].k
        nk = n - k
        H = np.zeros((nk, X_OR_Z * n), dtype="i1")
        for d in range(self._code_depth):
            if d not in self.H_depth:
                self.H_depth[d] = {}
            for i in range(self.code_instances[d].H.shape[0]):
                H[nind] = self.get_mother_operator(self.code_instances[d].H[i], d)
                self.H_depth[d][nind] = i
                nind += 1
        return H

    def BP_decode(self, syndrome):
        """
        BP_decoding:
        input: syndrome
        output: L^T
        """
        # Estimation for Logical error
        error_probability = self.bitwise_error_probability
        conc_syndrome = syndrome.copy()
        for d in reversed(range(self._code_depth)):
            beta1 = np.zeros(self.code_instances[d].nk, dtype="i1")
            for mother_ind, follower_ind in self.H_depth[d].items():
                beta1[follower_ind] = conc_syndrome[mother_ind]
            self.code_instances[d].set_channel_param(error_probability, iid=False)
            if not self.code_instances[d].BITWISE:
                self.code_instances[d].BITWISE = True
                print(
                    "Warning: Forced BITWISE = True in self.code_instances["
                    + str(d)
                    + "]."
                )
            c1_decoder_output = self.code_instances[d].decode(beta1, mode="ML")
            L, s, error_probability = (
                c1_decoder_output["L"],
                self.code_instances[d].get_syndrome({"E": c1_decoder_output["T"]}),
                c1_decoder_output["BIT_WISE_LOGICAL_ERROR_PROBABILITY"],
            )
            if d != 0:
                conc_syndrome ^= self.get_syndrome(
                    {"E": self.get_mother_operator(self.code_instances[d].get_T(s), d)}
                )
        L0 = self.get_mother_operator(L, 0)
        T = np.zeros(X_OR_Z * self.n, dtype="i1")
        for d in range(1, self._code_depth):
            T ^= self.get_mother_operator(self.code_instances[d].decoder_output["T"], d)
        self.decoder_output["L"] = L0
        self.decoder_output["T"] = self.get_T(syndrome)
        self.decoder_output["LT"] = L0 ^ T ^ self.get_T(conc_syndrome)
        self.decoder_output["LOGICAL_ERROR_PROBABILITY"] = error_probability

    def decode(self, syndrome, **param):
        # self._mode = "BP"
        if "BP" == self.mode:
            self.BP_decode(syndrome)
        else:
            super().decode(syndrome, **param)
        return self.decoder_output

    @property
    def code_instances(self):
        return self._code_instances

    @property
    def code_depth(self):
        return self._code_depth
