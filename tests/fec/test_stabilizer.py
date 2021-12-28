import numpy as np
import pytest
from src.fec import *

TEST_DATA = {
    "T_ind" : [[1,1,0,0],[1,0,0,1]],
    "T" : [[0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0]],
    "S_ind" : [3,0,15],
    "S" : [[1,1,1,1,0,1,0,0,1,0],[0,0,0,0,0,0,0,0,0,0],[0,0,1,0,1,1,1,0,0,0]],
    "L_ind" : [0,2,1,[0,0],[1,0],[0,1]],
    "L" : [[0,0,0,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0],[0,0,0,0,0,1,1,1,1,1],[0,0,0,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0],[0,0,0,0,0,1,1,1,1,1]],
}

def test_five():
    c = FIVE()
    for i in range(2 ** (c.n-c.k)):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i,4))))==i

    for t in range(len(TEST_DATA["T_ind"])):
        assert 0==sum(np.abs(c.get_T(TEST_DATA["T_ind"][t]) - np.array(TEST_DATA["T"][t])))

    for s in range(len(TEST_DATA["S_ind"])):
        assert 0==sum(np.abs(c.get_S(TEST_DATA["S_ind"][s]) - np.array(TEST_DATA["S"][s])))

    for l in range(len(TEST_DATA["L_ind"])):
        assert 0==sum(np.abs(c.get_L(TEST_DATA["L_ind"][l]) - np.array(TEST_DATA["L"][l])))

    for beta in range(2 ** c.k):
        assert 0==sum(np.abs(c.get_syndrome(c.get_L(beta))))

def test_STEANE():
    c = STEANE()
    for i in range(64):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i,6))))==i