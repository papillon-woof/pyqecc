import numpy as np
import pytest
from src.util import *
ep = 10 ** -10
def test_hard_decision_len():
    assert np.all(hard_decision_len(np.array([1.1,-1.2,3.1,4.1]))==np.array([0,1,0,0]).astype('i8'))
    assert np.all(hard_decision_len(np.array([0,-0,0,0]))==np.array([1,1,1,1]).astype('i8'))

def test_ideal_power():
    with pytest.raises(ValueError):
        ideal_power(0)
    s = 2*np.random.randint(0,4,100000)+1
    assert abs(ideal_power(1)-1)<ep
    assert abs(ideal_power(2)-5)<ep
    assert abs(ideal_power(3)-21)<ep
    assert abs(ideal_power(3)-power(s))<10 ** -1

def test_to_db():
    with pytest.raises(ValueError):
        to_db(-1)
    assert to_db(0) ==float("-inf")
    assert to_db(10) == 10
    assert to_db(1) == 0

def test_to_r():
    assert to_r(to_db(0)) == 0
    assert to_r(to_db(10)) == 10
    assert to_r(to_db(1)) == 1

def test_power():
    assert power(-1) == 1
    assert power([2,1]) == 2.5
    with pytest.raises(TypeError):
        power("1")
    assert power(np.array([2,1])) == 2.5

def test_awgn_channel():
    with pytest.raises(ValueError):
        awgn_channel(5)
    with pytest.raises(ValueError):
        awgn_channel(5,N0=-1)
    with pytest.raises(ValueError):
        awgn_channel(5,N=-5)
    assert sum(awgn_channel(5,N=0)) == 0
    assert abs(sum(awgn_channel(10000000,N=1) ** 2)/10000000-1) < 10 ** -2
    assert abs(sum(awgn_channel(10000000,N=10) ** 2)/10000000-10) < 10 ** -2
    assert abs(sum(awgn_channel(10000000,N0=5) ** 2)/10000000-2.5) < 10 ** -2
    assert abs(sum(awgn_channel(10000000,N0=12) ** 2)/10000000-6) < 10 ** -2
