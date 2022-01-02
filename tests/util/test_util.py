import numpy as np
import pytest
from src.util import *
DIFF_THRESHOLD = 1E-4
def test_symplex_binary_inner_product():
    assert 1==symplex_binary_inner_product(np.array([0,0,1,0],dtype='i1'),np.array([1,0,0,0],dtype='i1'))
    assert 0==symplex_binary_inner_product(np.array([1,0,0,0],dtype='i1'),np.array([1,0,0,0],dtype='i1'))
    assert 0==symplex_binary_inner_product(np.array([1,0,1,0],dtype='i1'),np.array([1,0,1,0],dtype='i1'))
    assert 1==symplex_binary_inner_product(np.array([1,1,1,0],dtype='i1'),np.array([1,0,0,0],dtype='i1'))

def test_arr2int():
    case = np.array([[0,0,0,0],[0,0,0,1],[0,0,1,0],[0,0,1,1],[0,1,0,0],[0,1,0,1],[0,1,1,0],[0,1,1,1],[1,0,0,0],[1,0,0,1],[1,0,1,0],[1,0,1,1],[1,1,0,0],[1,1,0,1],[1,1,1,0],[1,1,1,1]])
    for i in range(len(case)):
        assert i == arr2int(case[i])

def test_int2arr():
    m = 5
    for i in range(2**m):
        assert i== arr2int(int2arr(i,m))

TEST_DATA = {
    "QUBIT_NUMBER": [2,3,4,5,6,7]
}

def test_bitwise_to_blockwise_probability():
    for n in TEST_DATA["QUBIT_NUMBER"]:
        p = np.random.rand(4*n).reshape(n,4)
        for i in range(n):
            p[i]/=sum(p[i])
        assert DIFF_THRESHOLD>np.sum(np.abs(blockwise_to_bitwise_probability(bitwise_to_blockwise_probability(p))-p))