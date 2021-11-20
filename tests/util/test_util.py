import numpy as np
import pytest
from src.util import *
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