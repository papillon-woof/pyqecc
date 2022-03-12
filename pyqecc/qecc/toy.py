import numpy as np
from .stabilizer import *
from .concatenated import *


def SteaneCode(mode="ML_LUT"):
    N = 7
    K = 1
    H = np.array(
        [
            [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1],
        ],
        dtype="i1",
    )
    TX = {
        0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        1: np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        2: np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        3: np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        4: np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        5: np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        6: np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
        7: np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]),
    }
    TZ = {
        0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        1: np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]),
        2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),
        3: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]),
        4: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
        5: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]),
        6: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),
        7: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]),
    }
    T = {}
    for i in range(8):
        for j in range(8):
            T[8 * j + i] = TX[i] + TZ[j]
    L = np.array(
        [
            [
                [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
            ]
        ],
        dtype="i1",
    )
    sc = SC(N, K, H=H, T=T, L=L, mode=mode)
    sc._name = "STEANE_CODE"
    return sc


def FiveCode(mode="ML_LUT"):
    H = np.array(
        [
            [1, 0, 0, 1, 0, 0, 1, 1, 0, 0],
            [0, 1, 0, 0, 1, 0, 0, 1, 1, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 1, 1],
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 1],
        ],
        dtype="i1",
    )
    T = {
        0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        1: np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        8: np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
        12: np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0]),
        6: np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0]),
        3: np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),
        10: np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0]),
        5: np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
        2: np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0]),
        9: np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),
        4: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1]),
        11: np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0]),
        13: np.array([0, 1, 0, 0, 0, 0, 1, 0, 0, 0]),
        14: np.array([0, 0, 1, 0, 0, 0, 0, 1, 0, 0]),
        15: np.array([0, 0, 0, 1, 0, 0, 0, 0, 1, 0]),
        7: np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 1]),
    }
    L = np.array(
        [[[1, 1, 1, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]]], dtype="i1"
    )
    sc = SC(5, 1, H=H, T=T, L=L, mode=mode)
    sc._name = "FIVE_CODE"
    return sc


def BitFlipCode(mode="ML_LUT"):
    H = np.array([[0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 1, 1]], dtype="i1")
    L = np.array([[[1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1]]], dtype="i1")
    T = {
        0: np.array([0, 0, 0, 0, 0, 0]),
        1: np.array([0, 0, 1, 0, 0, 0]),
        2: np.array([1, 0, 0, 0, 0, 0]),
        3: np.array([0, 1, 0, 0, 0, 0]),
    }
    sc = SC(3, 1, H=H, T=T, L=L, mode=mode)
    sc._name = "BIT_FLIP_CODE"
    return sc


def PhaseFlipCode(mode="ML_LUT"):
    H = np.array([[1, 1, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0]], dtype="i1")
    T = {
        0: np.array([0, 0, 0, 0, 0, 0]),
        1: np.array([0, 0, 0, 0, 0, 1]),
        2: np.array([0, 0, 0, 1, 0, 0]),
        3: np.array([0, 0, 0, 0, 1, 0]),
    }
    L = np.array([[[1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1]]], dtype="i1")
    sc = SC(3, 1, H=H, T=T, L=L, mode=mode)
    sc._name = "PHASE_FLIP_FLIP_CODE"
    return sc


def ShorCode(mode="ML_LUT"):
    c0 = PhaseFlipCode()
    c1 = ParaCode([BitFlipCode for i in range(3)])
    sc = ConcCode([c0, c1])
    sc._name = "SHOR_CODE"
    return sc
