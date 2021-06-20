import numpy as np
from src import *

myst = STEANE()
print(myst)
e = np.zeros(14)

ee = myst.hard_decode(e)

def BLER(MONTE=10000):
    myst = STEANE()
    for p in [0.005*i for i in range(20)]:
        for monte in range(MONTE):
            e = depolarizing_noise(p)
            ee = myst.hard_decode(e)
