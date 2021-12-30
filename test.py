import numpy as np
from src import *
#c = CombCode([FIVE() for i in range(2)])
c1 = CombCode([BIT_FLIP(),BIT_FLIP(),BIT_FLIP()])
c0 = PHASE_FLIP()
c = ConcCode([c0,c1])
for i in range(2 ** 8):
    print(i,any2arr(i,8),c.get_T(any2arr(i,8)),c.get_syndrome(c.get_T(any2arr(i,8))))
exit()
dec_sim(c);exit()