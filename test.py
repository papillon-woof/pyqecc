import numpy as np
from src import *
c1 = CombCode([BIT_FLIP(),BIT_FLIP(),BIT_FLIP()])
c0 = PHASE_FLIP()
c = ConcCode([c0,c1])
dec_sim(c1);exit()
p = 0.1
c.set_block_wise_p(np.array([1-p,p/3,p/3,p/3]))
s = np.zeros(8)
print(c.BP_decode(s))