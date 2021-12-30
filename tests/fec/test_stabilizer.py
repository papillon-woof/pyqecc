import numpy as np
import pytest
from src.fec import *
DIFF_THRESHOLD = 1E-4

TEST_FIVE_DATA = {
    "T_ind" : [[1,1,0,0],[1,0,0,1]],
    "T" : [[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0]],
    "S_ind" : [3,0,15],
    "S" : [[1,1,1,1,0,1,0,0,1,0],[0,0,0,0,0,0,0,0,0,0],[0,0,1,0,1,1,1,0,0,0]],
    "L_ind" : [0,2,1,[0,0],[1,0],[0,1]],
    "L" : [[0,0,0,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0],[0,0,0,0,0,1,1,1,1,1],[0,0,0,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0],[0,0,0,0,0,1,1,1,1,1]],
}

def test_general():
    for c in [FIVE(),STEANE(),BIT_FLIP(),PHASE_FLIP()]:
        # test for Error
        for e in range(2 ** (X_OR_Z*c.n)):
            ee = any2arr(e,X_OR_Z*c.n)
            s = c.get_syndrome(ee)
            assert sum(np.abs(s-c.get_syndrome(c.get_T(s))))==0

def test_set_error_probability():
    c = BIT_FLIP()
    P_BITWISE = np.array([[1,0,0,0],[0.9,0.2,0.1,0],[1,0,0,0]])
    P_IID = np.array([0.7,0.1,0.1,0.1])
    P_BLOCKWISE = np.random.rand(2 ** (2*3));P_BLOCKWISE/=sum(P_BLOCKWISE)
    print(P_BLOCKWISE)
    c.set_error_probability(P_BITWISE,iid=False,BITWISE=True)
    assert np.abs(c.get_error_probability([0,0,0,0,0,0])-0.9)<DIFF_THRESHOLD
    assert np.abs(c.get_error_probability([0,0,1,0,0,1])-0)<DIFF_THRESHOLD
    assert np.abs(c.get_error_probability([0,1,0,0,0,0])-0.1)<DIFF_THRESHOLD
    c.set_error_probability(P_IID,BITWISE=True)
    assert np.abs(c.get_error_probability([0,0,0,0,0,0])-(0.7 ** 3))<DIFF_THRESHOLD
    assert np.abs(c.get_error_probability([0,0,1,0,0,1])-(0.7 ** 2)*0.1)<DIFF_THRESHOLD
    assert np.abs(c.get_error_probability([0,1,0,0,0,0])-(0.7 ** 2)*0.1)<DIFF_THRESHOLD
    c.set_error_probability(P_BLOCKWISE,BITWISE=False)
    assert np.abs(c.get_error_probability([0,0,0,0,0,0])-P_BLOCKWISE[arr2int([0,0,0,0,0,0])])<DIFF_THRESHOLD
    assert np.abs(c.get_error_probability([0,0,1,0,0,1])-P_BLOCKWISE[arr2int([0,0,1,0,0,1])])<DIFF_THRESHOLD
    assert np.abs(c.get_error_probability([0,1,0,0,0,0])-P_BLOCKWISE[arr2int([0,1,0,0,0,0])])<DIFF_THRESHOLD

def test_five():
    c = FIVE()
    for i in range(2 ** (c.n-c.k)):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i,4))))==i

    for t in range(len(TEST_FIVE_DATA["T_ind"])):
        assert 0==sum(np.abs(c.get_T(TEST_FIVE_DATA["T_ind"][t]) - np.array(TEST_FIVE_DATA["T"][t])))

    for s in range(len(TEST_FIVE_DATA["S_ind"])):
        assert 0==sum(np.abs(c.get_S(TEST_FIVE_DATA["S_ind"][s]) - np.array(TEST_FIVE_DATA["S"][s])))

    for l in range(len(TEST_FIVE_DATA["L_ind"])):
        assert 0==sum(np.abs(c.get_L(TEST_FIVE_DATA["L_ind"][l]) - np.array(TEST_FIVE_DATA["L"][l])))

    for beta in range(2 ** c.k):
        assert 0==sum(np.abs(c.get_syndrome(c.get_L(beta))))

def test_STEANE():
    c = STEANE()
    for i in range(64):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i,6))))==i