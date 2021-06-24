import numpy as np
import pytest
from src.fec import *

def test_toy():
    c = FIVE()
    assert 0==sum(c.get_T([1,1,0,0]) - np.array([1,1,0,0,0,0,0,0,0,0]))
    assert 0==sum(c.get_T([1,0,0,1]) - np.array([1,0,0,0,0,0,0,0,0,1]))

    assert 0==sum(c.get_S(3) - np.array([1,1,1,1,0,1,0,0,0,1]))
    assert 0==sum(c.get_S(0) - np.array([0,0,0,0,0,0,0,0,0,0]))

    assert 0==sum(c.get_L(0) - np.array([0,0,0,0,0,0,0,0,0,0]))
    assert 0==sum(c.get_L(1) - np.array([1,1,1,1,1,0,0,0,0,0]))
    assert 0==sum(c.get_L(2) - np.array([0,0,0,0,0,1,1,1,1,1]))
    assert 0==sum(c.get_syndrome(c.get_L(0)) - np.zeros(4,dtype='i1'))
    assert 0==sum(c.get_syndrome(c.get_L(1)) - np.zeros(4,dtype='i1'))
    assert 0==sum(c.get_syndrome(c.get_L(2)) - np.zeros(4,dtype='i1'))
    assert 0==sum(c.get_syndrome(c.get_L(3)) - np.zeros(4,dtype='i1'))
