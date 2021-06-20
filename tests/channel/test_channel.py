import numpy as np
import pytest
from src.channel import *
def test_depolarizing():
    np.random.seed=1
    n = 100000
    for p in [0.01,0.1,0.2]:
        e = depolarizing_noise(n,p)
        cx = 0;cy = 0;cz = 0
        for i in range(n):
            if e[i]==1 and e[i+n]==1:
                cy += 1
            elif e[i]==1:
                cx += 1
            elif e[i+n]==1:
                cz += 1
        #print(cy,cx,cz)
        #print(cy/n,cx/n,cz/n)
        assert np.abs(cy/n-p/3)<0.003
        assert np.abs(cx/n-p/3)<0.003
        assert np.abs(cz/n-p/3)<0.003
