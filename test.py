import numpy as np
from src import *
c2 = ParaCode([FIVE() for i in range(25)])
c1 = ParaCode([FIVE() for i in range(5)])
c0 = FIVE()
ccc = ConcCode([c0,c1,c2])
cc = ConcCode([c0,c1])
e = np.zeros(2*cc.n)
np.set_printoptions(threshold=100000)
prob = [0.13,0.15,0.1885]
dec_sim(c0,PROB=prob,MONTE=1000,LOG_OUTPUT_SPAN=10)
dec_sim(cc,PROB=prob,MONTE=1000,LOG_OUTPUT_SPAN=10)
dec_sim(ccc,PROB=prob,MONTE=1000,LOG_OUTPUT_SPAN=10);exit()
