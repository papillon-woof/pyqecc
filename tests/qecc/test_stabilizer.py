import numpy as np
import pytest
from pyqecc.qecc import *

DIFF_THRESHOLD = 1e-4

TEST_FIVE_DATA = {
    "T_ind": [[1, 1, 0, 0], [1, 0, 0, 1]],
    "T": [[0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]],
    "S_ind": [3, 0, 15],
    "S": [
        [1, 1, 1, 1, 0, 1, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 1, 1, 1, 0, 0, 0],
    ],
    "L_ind": [0, 2, 1, [0, 0], [1, 0], [0, 1]],
    "L": [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    ],
}


def test_general():
    for c in [FiveCode(), SteaneCode(), BitFlipCode(), PhaseFlipCode()]:
        # test for Error
        for e in range(2 ** (X_OR_Z * c.n)):
            ee = any2arr(e, X_OR_Z * c.n)
            s = c.get_syndrome(ee)
            assert sum(np.abs(s - c.get_syndrome(c.get_T(s)))) == 0


TEST_SET_ERROR_PROBABILITY_DATA = {
    "NUM_OF_CASE": 3,
    "P_BITWISE": np.array([[1, 0, 0, 0], [0.9, 0.2, 0.1, 0], [1, 0, 0, 0]]),
    "TEST_CASE_P_BITWISE": [[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 1, 0, 0, 0, 0]],
    "RESULT_P_BITWISE": [0.9, 0, 0.1],
    "P_IID": np.array([0.7, 0.1, 0.1, 0.1]),
    "TEST_CASE_P_IID": [[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 1, 0, 0, 0, 0]],
    "RESULT_P_IID": [(0.7 ** 3), (0.7 ** 2) * 0.1, (0.7 ** 2) * 0.1],
    "P_BLOCKWISE": np.random.rand(2 ** (2 * 3)),
    "TEST_CASE_P_BLOCKWISE": [
        [0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1],
        [0, 1, 0, 0, 0, 0],
    ],
    # "RESULT_P_BLOCKWISE" : [0.9,0,0.1],
}


def test_set_error_probability():
    c = BitFlipCode()
    c.set_error_probability(
        TEST_SET_ERROR_PROBABILITY_DATA["P_BITWISE"], iid=False, BITWISE=True
    )
    for i in range(TEST_SET_ERROR_PROBABILITY_DATA["NUM_OF_CASE"]):
        assert (
            np.abs(
                c.get_error_probability(
                    TEST_SET_ERROR_PROBABILITY_DATA["TEST_CASE_P_BITWISE"][i]
                )
                - TEST_SET_ERROR_PROBABILITY_DATA["RESULT_P_BITWISE"][i]
            )
            < DIFF_THRESHOLD
        )
    c.set_error_probability(TEST_SET_ERROR_PROBABILITY_DATA["P_IID"], BITWISE=True)
    for i in range(TEST_SET_ERROR_PROBABILITY_DATA["NUM_OF_CASE"]):
        assert (
            np.abs(
                c.get_error_probability(
                    TEST_SET_ERROR_PROBABILITY_DATA["TEST_CASE_P_IID"][i]
                )
                - TEST_SET_ERROR_PROBABILITY_DATA["RESULT_P_IID"][i]
            )
            < DIFF_THRESHOLD
        )

    TEST_SET_ERROR_PROBABILITY_DATA["P_BLOCKWISE"] /= sum(
        TEST_SET_ERROR_PROBABILITY_DATA["P_BLOCKWISE"]
    )
    c.set_error_probability(
        TEST_SET_ERROR_PROBABILITY_DATA["P_BLOCKWISE"], BITWISE=False
    )
    for i in range(TEST_SET_ERROR_PROBABILITY_DATA["NUM_OF_CASE"]):
        ind = TEST_SET_ERROR_PROBABILITY_DATA["TEST_CASE_P_BLOCKWISE"][i]
        assert (
            np.abs(
                c.get_error_probability(ind)
                - TEST_SET_ERROR_PROBABILITY_DATA["P_BLOCKWISE"][arr2int(ind)]
            )
            < DIFF_THRESHOLD
        )


def test_FiveCode():
    c = FiveCode()
    for i in range(2 ** (c.n - c.k)):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i, 4)))) == i
    for t in range(len(TEST_FIVE_DATA["T_ind"])):
        assert 0 == sum(
            np.abs(
                c.get_T(TEST_FIVE_DATA["T_ind"][t]) - np.array(TEST_FIVE_DATA["T"][t])
            )
        )
    for s in range(len(TEST_FIVE_DATA["S_ind"])):
        assert 0 == sum(
            np.abs(
                c.get_S(TEST_FIVE_DATA["S_ind"][s]) - np.array(TEST_FIVE_DATA["S"][s])
            )
        )
    for l in range(len(TEST_FIVE_DATA["L_ind"])):
        assert 0 == sum(
            np.abs(
                c.get_L(TEST_FIVE_DATA["L_ind"][l]) - np.array(TEST_FIVE_DATA["L"][l])
            )
        )
    for beta in range(2 ** c.k):
        assert 0 == sum(np.abs(c.get_syndrome(c.get_L(beta))))


def test_SteaneCode():
    c = SteaneCode()
    for i in range(64):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i, 6)))) == i
