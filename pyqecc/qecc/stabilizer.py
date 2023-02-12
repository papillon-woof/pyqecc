from typing import List, Dict, Union
import warnings
import numpy as np
import numpy.typing as npt
from .abstruct import CODE
from ..util import (
    bitwise_to_blockwise_error_probability,
    blockwise_to_bitwise_error_probability,
    arr2int,
    symplex_binary_inner_product,
    any2arr,
    X_OR_Z,
    gaussjordan,
    int2arr,
)


class SC(CODE):
    """
    Class of stabilizer code

    Attributes
    ----------
    n: int
        Codeword length
    k: int
        Information length
    H: numpy.ndarray[numpy.int8]
        Parity check matrix
    T: numpy.ndarray[numpy.int8]
        Recovery operator
    L: numpy.ndarray[numpy.int8]
        Logical operator
    mode: str
        Mode of decoder
    BITWISE: bool
        Calculation by using bitwise or blockwise
    ANALOG_INFORMATON: str
        Selector of using Analog information
    """

    NAME = "stabilizer code"
    USAGE = 'The argments for class of stabilizer code\'s are the codeword length "n" and information length "k".'

    def __init__(
        self,
        n: int,
        k: int,
        H: npt.NDArray[np.int8],
        T: npt.NDArray[np.int8],
        L: npt.NDArray[np.int8],
        P=None,
        mode: str="HD",
        BITWISE: bool=True,
        ANALOG_INFORMATION: str="No",
    ) -> None:
        """Constractor
        
        Parameters
        ----------
        n: int
            Codeword length
        k: int
            Information length
        _H: numpy.ndarray[numpy.int8]
            Parity check matrix
        _T: numpy.ndarray[numpy.int8]
            Recovery operator
        _L: numpy.ndarray[numpy.int8]
            Logical operator
        _LUT: 
            Lookup-table for multiplication result of recovery operator and logical operator
        mode: str
            Mode of decoder
        BITWISE: bool
            Selecter whether using bitwise or blockwise
        ANALOG_INFORMATON: str
            Selecter of using Analog information
        """
        self.decoder_output = {
            "LT": None,
            "L": None,
            "T": None,
            "LOGICAL_ERROR_PROBABILITY": None,
            "BIT_WISE_LOGICAL_ERROR_PROBABILITY": None,
        }
        super().__init__(n, k)
        self._ML_DECODING_QUBITS_LIMIT = 10
        self._mode = mode
        self._H = H
        self._T = T
        self._L = L
        self._BITWISE = BITWISE
        self.set_channel_param(P, self.BITWISE)
        self._LUT: Dict[int,npt.NDArray[np.int8]] = {}
        self._blockwise_error_probability = False
        self._bitwise_error_probability = False

    def set_channel_param(
        self, error_probability:npt.NDArray[np.float64], BITWISE: bool=True, iid: bool=True, warnings_ignore: bool=False
    )  -> None:
        """
        Parametors
        ----------
        error_probability: numpy.ndarray(numpy.float64)
            Channel side information
        BITWISE:
            Selecter whether using bitwise or blockwise
        iid:
            Assuming the i.i.d (Independent and identically distributed) channel
        warnings_ignore: bool
            Selecter whether output warning-message or not
        """
        if not error_probability is None:
            if BITWISE:
                if iid:
                    error_probability = np.array(
                        [error_probability for i in range(self.n)]
                    )
                self.bitwise_error_probability = error_probability
                if self.ML_DECODING_QUBITS_LIMIT < self.n:
                    if warnings_ignore:
                        warnings.warn(
                            "Warning: The num. of qubit exceed ML_DECODING_QUBITS_LIMIT."
                        )
                    self.blockwise_error_probability = False
                else:
                    self.blockwise_error_probability = (
                        bitwise_to_blockwise_error_probability(error_probability)
                    )  # By approximate the bitwise probability.
            else:
                self.blockwise_error_probability = error_probability
                self.bitwise_error_probability = blockwise_to_bitwise_error_probability(
                    error_probability
                )
        else:
            if warnings_ignore:
                warnings.warn("Warning: The input error probability doesn't set to variable.")

    def get_error_probability(self, E: npt.NDArray[np.int8]):
        """
        Getter of error probability (channel side information)

        Parameters
        ----------
        E: numpy.ndarray(numpy.int8)
            error operator
        """
        return self.blockwise_error_probability[arr2int(E)]

    def get_syndrome(self, channel_output: Dict):
        """
        Getter of syndrome

        Parameters
        ----------
        channel_output: dict
            channel 
        """

        return symplex_binary_inner_product(self._H, channel_output["E"])

    def get_T(self, ind: int):
        """
        Getter of recovery operator

        Parameters
        ----------
        ind: dict
            channel 
        """
        return self.T[arr2int(ind)].astype("i1")  # LUTでの計算? BPでの計算もあり?LDPCについて学ぶ．

    def get_S(self, ind_list: List[int]):
        S = np.zeros(X_OR_Z * self.n, dtype="i1")
        ind_list = any2arr(ind_list, self.nk)
        for i in range(self.nk):
            S ^= ind_list[i] & self.H[i]
        return S

    def get_L(self, alpha: npt.NDArray[np.int8]):
        """

        parameters
        ----------
        alpha: numpy.ndarray of numpy.int8
            Vector_representation of operator
        
        returns
        -------
        L: numpy.ndarray of numpy.int8
            Vector representation of logical operator
            [LX1,LX2,LX3,...,LX(n-k)|LZ1,LZ2,LZ3,...,LZ(n-k)]
        """
        alpha = any2arr(alpha, X_OR_Z * self.k)
        L = np.zeros(X_OR_Z * self.n, dtype="i1")
        # LX or LZ
        for j in range(self.k):
            for i in range(X_OR_Z):
                L ^= alpha[i * self.k + j] & self.L[j][i]
        return L

    def in_S(self, b: npt.NDArray[np.int8]) -> bool:
        """
        Return whether S include b or not

        parameters
        ----------
        b : Vector representation
        """
        return (
            sum(
                gaussjordan(np.append(self.H, b).reshape(self.nk + 1, X_OR_Z * self.n))[
                    self.nk
                ]
            )
            == 0
        )

    def set_LUT(self) -> None:
        for i in range(2 ** (self.nk)):
            self._LUT[i] = self.ML_decode(int2arr(i, (self.nk)))

    def ML_decode(self, syndrome, **param):
        """
        Maximum likelihood decoding.

        Parameters
        ----------
        syndrome: numpy.ndarray of numpy.int8
            Syndrome
        """
        if self.n > self.ML_DECODING_QUBITS_LIMIT:
            raise ValueError(
                "Error: The qubit n ="
                + str(self.n)
                + " is limited because of a large decoding complexity. You can change the qubit limit."
            )
        T = self.get_T(syndrome)
        logical_error_probability = np.zeros(2 ** (X_OR_Z * self.k))
        for lind in range(2 ** (X_OR_Z * self.k)):
            L = self.get_L(lind)
            for sind in range(2 ** (self.nk)):
                S = self.get_S(sind)
                E = L ^ T ^ S
                logical_error_probability[lind] += self.get_error_probability(E)
        llind = np.argmax(logical_error_probability)
        self.decoder_output["L"] = self.get_L(llind)
        self.decoder_output["T"] = self.get_T(syndrome)
        self.decoder_output["LT"] = self.decoder_output["L"] ^ self.decoder_output["T"]
        self.decoder_output["LOGICAL_ERROR_PROBABILITY"] = logical_error_probability
        return self.decoder_output["LT"]

    def hard_decode(self, syndrome) -> None:
        """
        Hard decision decoding (HDD)
        Parameters
        ----------
        syndrome: numpy.ndarray of numpy.int8
            Syndrome
        """
        self.decoder_output["LT"] = self.get_T(syndrome)

    def LUT_decode(self, syndrome) -> None:
        """
        ML decoding by using LUT [Delete or fix in the future]

        Parameters
        ----------
        syndrome: numpy.ndarray of numpy.int8
            Syndrome
        """
        self.decoder_output["LT"] = self.LUT[arr2int(syndrome)]

    def decode(
        self,
        syndrome: npt.NDArray[np.int8],
        mode: str="",
        **param: Dict,
    ) -> Dict:
        """
        decoder

        Parameters
        ----------
        syndrome: numpy.ndarray of numpy.int8
            Syndrome
        mode: str
            Mode of decoder
        **param: dict
            Optional parameter 
        """
        if not mode is False:
            self._mode = mode
        if self._mode == "ML":
            self.ML_decode(syndrome, **param)
        if self._mode == "ML_LUT":
            if self.LUT == {}:
                self.set_LUT()
            self.LUT_decode(syndrome, **param)
        if self._mode == "HD":
            self.hard_decode(syndrome, **param)
        if self._mode == "ANALOG":
            self.analog_decode(syndrome[0], syndrome[1])
        return self.decoder_output

    @property
    def blockwise_error_probability(self):
        return self._blockwise_error_probability

    @blockwise_error_probability.setter
    def blockwise_error_probability(self, blockwise_error_probability):
        self._blockwise_error_probability = blockwise_error_probability

    @property
    def bitwise_error_probability(self):
        return self._bitwise_error_probability

    @bitwise_error_probability.setter
    def bitwise_error_probability(self, bitwise_error_probability):
        self._bitwise_error_probability = bitwise_error_probability

    @property
    def L(self) -> npt.NDArray[np.int8]:
        return self._L

    @property
    def T(self) -> npt.NDArray[np.int8]:
        return self._T

    @property
    def H(self) -> npt.NDArray[np.int8]:
        return self._H

    @property
    def LUT(self) -> Dict[int,npt.NDArray[np.int8]]:
        return self._LUT

    @property
    def mode(self) -> str:
        return self._mode

    @property
    def ML_DECODING_QUBITS_LIMIT(self) -> int:
        return self._ML_DECODING_QUBITS_LIMIT

    @property
    def BITWISE(self) -> bool:
        return self._BITWISE

    @BITWISE.setter
    def BITWISE(self, BITWISE):
        self._BITWISE = BITWISE

    def __str__(self):
        output = ""
        output += "NAME               : " + str(self.name) + "\n"
        output += "PHYSICAL QUBITS (n): " + str(self.n) + "\n"
        output += "LOGICAL QUBITS (k) : " + str(self.k) + "\n"
        output += "CODE RATE (R = k/n): " + str(self.R) + "\n"
        output += "DECODING_MODE      : " + str(self.mode) + "\n"
        return output
