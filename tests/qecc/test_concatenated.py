import numpy as np
import pytest
from pyqecc.qecc import *

TEST_DATA_PARA = [
    {
    "TEST_NUM": 3,
    "PARA_NUM": 1,
    "T": [[0,0,0,0],[1,1,0,0],[1,0,0,1]],
    "L": [[0,0],[1,0],[0,1]],
    "GET_T" : [np.array([0,0,0,0,0,0,0,0,0,0]),np.array([0,0,1,0,0,0,0,0,0,0]),np.array([0,0,0,0,0,0,0,0,1,0])],
    "GET_L" : [np.array([0,0,0,0,0,0,0,0,0,0]),np.array([1,1,1,1,1,0,0,0,0,0]),np.array([0,0,0,0,0,1,1,1,1,1])],
    },
    {
    "TEST_NUM": 3,
    "PARA_NUM": 2,
    "T": [[0,0,0,0,0,0,0,0],[1,1,0,0,0,0,0,0],[1,1,0,0,1,0,0,1]],
    "L": [[0,0,0,0],[1,0,0,0],[0,1,1,0]],
    "GET_T" : [np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])],
    "GET_L" : [np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),np.array([1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),np.array([0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0])],
    },
]

TEST_DATA_CONC = [
    {
    "TEST_NUM_T": 4,
    "TEST_NUM_L": 6,
    "TEST_CODE": ConcCode([PhaseFlipCode(),ParaCode([BitFlipCode(),BitFlipCode(),BitFlipCode()])]),
    "T": [[0,0,0,0,0,0,0,0],[1,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,1],[1,0,0,0,1,0,0,1]],
    "L": [[0,0],[1,0],[0,1],[1,1],2,1],
    "GET_T" : [np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),np.array([0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0]),np.array([0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]),np.array([0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0]),np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),np.array([0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0])],
    "GET_L" : [np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),np.array([1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]),np.array([0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1]),np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),np.array([1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]),np.array([0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1])],
    },
    {
    "TEST_NUM_T": 1,
    "TEST_NUM_L": 6,
    "TEST_CODE": ConcCode([FiveCode(),ParaCode([FiveCode() for i in range(5)])]),
    "T": [
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        ],
    "L": [[0,0],[1,0],[0,1],[1,1],2,1],
    "GET_T" : [
        np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
        ],
    "GET_L" : [
        np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
        np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
        np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),
        np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),
        np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
        np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),
        ],
    },
]

TEST_DATA_CONC_H = [
    {
    "TEST_NUM_S": 6,
    "TEST_CODE": [
        ConcCode([FiveCode(),ParaCode([FiveCode() for i in range(5)])]),
        ConcCode([FiveCode(),ParaCode([FiveCode() for i in range(5)]),ParaCode([FiveCode() for i in range(25)])]),
        ],
    "GET_S" : [[2,4],[1],[3,5],[9,15],[1,3,5],[6,9,15],
        ],
    }
]

def test_CONB():
    for test_data in TEST_DATA_PARA:
        c = ParaCode([FiveCode() for i in range(test_data["PARA_NUM"])])
        for i in range(test_data["TEST_NUM"]):
            assert 0==sum(np.abs(c.get_T(test_data["T"][i]) - test_data["GET_T"][i]))
            assert 0==sum(np.abs(c.get_L(test_data["L"][i]) - test_data["GET_L"][i]))

def test_CONC():
    for test_data in TEST_DATA_CONC:
        c = test_data["TEST_CODE"]
        for i in range(test_data["TEST_NUM_T"]):
            assert 0==sum(np.abs(c.get_T(test_data["T"][i]) - test_data["GET_T"][i]))
        for i in range(test_data["TEST_NUM_L"]):
            assert 0==sum(np.abs(c.get_L(test_data["L"][i]) - test_data["GET_L"][i]))

def test_CONC_H():
    for test_data in TEST_DATA_CONC_H:
        for c in test_data["TEST_CODE"]:
            e = np.zeros(2*c.n)
            for j in test_data["GET_S"]:
                e[j]=1
            assert not c.in_S(e)
