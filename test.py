import numpy as np
from src import *
print(1&1);exit()
c2 = ParaCode([FIVE() for i in range(25)])
c1 = ParaCode([FIVE() for i in range(5)])
c0 = FIVE()
cc = ConcCode([c0,c1,c2])
ccc = ConcCode([c0,c1])
e = np.zeros(2*cc.n)
np.set_printoptions(threshold=100000)
dec_sim(ccc,PROB=[0.1885],MONTE=1000,LOG_OUTPUT_SPAN=10);exit()
