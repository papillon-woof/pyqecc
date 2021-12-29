import numpy as np
from src import *
c1 = CombCode([BIT_FLIP(),BIT_FLIP(),BIT_FLIP()])
c0 = PHASE_FLIP()
c = ConcCode([c0,c1])
print(c.get_T([1,1,0,0,0,0,0,0]))
print(c.get_T([1,0,0,0,1,0,0,1]))
