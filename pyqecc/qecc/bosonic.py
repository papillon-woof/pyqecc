import warnings
import numpy as np
from .abstruct import *
from ..util import *
from .stabilizer import SC


class GKP(CODE):
    NAME = "GKP qubit"
    USAGE = "GKP code class. This class uses the analog information (See ref )"
    START_MAX = -1e300

    def __init__(
        self,
        code_instance,
        sigma=0.1,
        mode="SYNDROME",
    ):
        self.decoder_output = {
            "LT": None,
            "LOGICAL_ERROR_PROBABILITY": None,
        }
        self._mode = mode
        self._code_instance = code_instance
        self._sigma = sigma
        super().__init__(self.code_instance.n, self.code_instance.k)
        # self.matrix_for_genarating_LLR = util.gaussjordan(np.concatenate([self.code_instance.H,self.code_instance.L],0),change=True)

    def set_channel_param(self, sigma):
        self._sigma = sigma

    def get_syndrome(self, channel_output):
        """
        引数: E: アナログ雑音もしくは誤り
        Eは，実際に生じたシフトである．
        """
        err = np.zeros(2 * self.n, dtype="i1")
        if "DELTA" not in channel_output:
            if channel_output["E"].dtype != "i1":
                raise ValueError("E is not binary. Please Chack the channel output")
            else:
                warnings.warn(
                    "Not include the analog information. Instedly, the decoder performs only discrete error because channel output doesn't include the analog information."
                )
                return self.code_instance.get_syndrome(channel_output)

        # 2√π>|E|>√π => error if 2√π < or < -2√π => +-2√π
        delta = pishifts(channel_output["DELTA"])
        e_pos = np.where(np.abs(delta) >= np.sqrt(np.pi) / 2)[0]
        err[e_pos] = 1
        syndrome = self.code_instance.get_syndrome({"E": err})

        # calculation Δm
        delta_m = np.abs(delta.copy())
        delta_m[e_pos] = (np.sqrt(np.pi) - delta_m)[e_pos]
        return syndrome, delta_m

    def in_S(self, s):
        return self.code_instance.in_S(s)

    def calc_llr(self, val, sigma):
        return np.log(
            np.exp(-(val**2) / (2 * sigma**2))
            / np.exp(-((np.sqrt(np.pi) - np.abs(val)) ** 2) / (2 * sigma**2))
        )

    def analog_ML_decode(self, information):
        """
        Calculate the maximum likelihood decoding (ML-decoding) which maximize the P(L|information)
        """
        syndrome = information[0]
        delta_m = information[1]
        most_likely_error = np.zeros(self.n, dtype="i1")
        ma = self.START_MAX

        for i in range(2 * 2**self.k):
            lt = self.code_instance.get_T(syndrome)
            for ii in range(2 * self.k):
                lt ^= (1 & (i >> ii)) * self.code_instance.get_L(i)
            for j in range(2 * 2 ** (self.n - self.k)):
                c = lt.copy()
                for jj in range(2 * self.nk):
                    c ^= (1 & (j >> jj)) * self.code_instance.get_S(j)
                llr = 0
                for k in range(len(c)):
                    llr += ((-1) ** c[k]) * self.calc_llr(delta_m[k], self.sigma)
                if ma < llr:
                    ma = llr
                    most_likely_error = c
        self.decoder_output["LT"] = most_likely_error
        return self.decoder_output

    def analog_decode(self, information):
        """
        Approximatly calculate the maximum likelihood decoding (ML-decoding)
        by the P(L|information) without degenerate for error.
        That is, stabilizer is not considering in calculation.
        """
        syndrome = information[0]
        delta_m = information[1]
        error = np.zeros(self.n, dtype="i1")
        ma = self.START_MAX
        for i in range(2 * 2**self.k):

            # About recovery opelator T
            lt = self.code_instance.get_T(syndrome)

            # About logical opelator L
            for ii in range(2 * self.k):
                lt ^= (1 & (i >> ii)) * self.code_instance.get_L(i)

            # Calculation for log-likelihood ratio (LLR)
            llr = 0
            for k in range(len(lt)):
                llr += ((-1) ** lt[k]) * self.calc_llr(delta_m[k], self.sigma)
            if ma < llr:
                ma = llr
                error = lt
        self.decoder_output["LT"] = error
        return self.decoder_output

    def decode(self, syndrome, mode=None):
        if mode is not None:
            self._mode = mode

        # Do not use the analog information
        if len(syndrome) == 1 or self.mode[:8] == "DIGITAL_":
            self.code_instance.decode(syndrome, mode=self.mode[8:])
            return super().decode(syndrome, analog=False)

        # Use the analog information
        elif self.mode == "ML":
            return self.analog_ML_decode(syndrome)
        elif self.mode == "SYNDROME":
            return self.analog_decode(syndrome)
        else:
            raise ValueError("Chack the decoder mode")

    @property
    def code_instance(self):
        return self._code_instance

    @property
    def sigma(self):
        return self._sigma

    @property
    def mode(self):
        return self._mode

    def __str__(self):
        output = super().__str__()
        output += "MODE            :" + str(self.mode) + "\n"
        return output
