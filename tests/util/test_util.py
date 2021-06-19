import numpy as np
import pytest
from src.util import *
def test_symplex_binary_inner_product():
    assert 1==symplex_binary_inner_product(np.array([0,0,1,0],dtype='i1'),np.array([1,0,0,0],dtype='i1'))
    assert 0==symplex_binary_inner_product(np.array([1,0,0,0],dtype='i1'),np.array([1,0,0,0],dtype='i1'))
    assert 0==symplex_binary_inner_product(np.array([1,0,1,0],dtype='i1'),np.array([1,0,1,0],dtype='i1'))
    assert 1==symplex_binary_inner_product(np.array([1,1,1,0],dtype='i1'),np.array([1,0,0,0],dtype='i1'))
    assert 0==symplex_binary_inner_product(np.array([[1,0,1,0],[0,0,1,1],[0,1,1,0]],dtype='i1'),np.array([1,0,1,0],dtype='i1'))[0]
    assert 1==symplex_binary_inner_product(np.array([[1,0,1,0],[0,0,1,1],[0,1,1,0]],dtype='i1'),np.array([1,0,1,0],dtype='i1'))[1]
    assert 1==symplex_binary_inner_product(np.array([[1,0,1,0],[0,0,1,1],[0,1,1,0]],dtype='i1'),np.array([1,0,1,0],dtype='i1'))[2]
