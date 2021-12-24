import numpy as np
import pytest
from src.fec import *

def test_toy():
    c = FIVE()
    assert 0==sum(np.abs(c.get_T([1,1,0,0]) - np.array([0,0,1,0,0,0,0,0,0,0])))
    assert 0==sum(np.abs(c.get_T([1,0,0,1]) - np.array([0,0,0,0,0,0,0,0,1,0])))

    assert 0==sum(np.abs(c.get_S(3) - np.array([1,1,1,1,0,1,0,0,1,0])))
    assert 0==sum(np.abs(c.get_S(0) - np.array([0,0,0,0,0,0,0,0,0,0])))
    assert 0==sum(np.abs(c.get_S(15) - np.mod(np.array([1,0,0,1,0,0,1,1,0,0])+np.array([0,1,0,0,1,0,0,1,1,0])+np.array([1,0,1,0,0,0,0,0,1,1])+np.array([0,1,0,1,0,1,0,0,0,1]),2)))
    
    assert 0==sum(np.abs(c.get_L(0) - np.array([0,0,0,0,0,0,0,0,0,0])))
    assert 0==sum(np.abs(c.get_L(2) - np.array([1,1,1,1,1,0,0,0,0,0])))
    assert 0==sum(np.abs(c.get_L(1) - np.array([0,0,0,0,0,1,1,1,1,1])))
    assert 0==sum(np.abs(c.get_L([0,0]) - np.array([0,0,0,0,0,0,0,0,0,0])))
    assert 0==sum(np.abs(c.get_L([1,0]) - np.array([1,1,1,1,1,0,0,0,0,0])))
    assert 0==sum(np.abs(c.get_L([0,1]) - np.array([0,0,0,0,0,1,1,1,1,1])))
    assert 0==sum(np.abs(c.get_syndrome(c.get_L(0)) - np.zeros(4,dtype='i1')))
    assert 0==sum(np.abs(c.get_syndrome(c.get_L(1)) - np.zeros(4,dtype='i1')))
    assert 0==sum(np.abs(c.get_syndrome(c.get_L(2)) - np.zeros(4,dtype='i1')))
    assert 0==sum(np.abs(c.get_syndrome(c.get_L(3)) - np.zeros(4,dtype='i1')))

def test_five():
    c = FIVE()
    for i in range(15):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i,4))))==i

def test_STEANE():
    c = STEANE()
    for i in range(64):
        assert arr2int(c.get_syndrome(c.get_T(int2arr(i,6))))==i