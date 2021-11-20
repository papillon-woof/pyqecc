import numpy as np
from src import *
c = STEANE()
for i in range(2 ** 6):
    print(c.get_T(int2arr(i,6)))

n = 100000
print(BLER(STEANE(mode='ML'),monte=10000))
print(BLER(STEANE(mode='ML'),monte=10000))
