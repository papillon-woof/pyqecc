import numpy as np
from src import *
c1 = CombCode([FIVE() for i in range(5)])
c0 = FIVE(mode='ML')
print(c0.get_syndrome([1,1,0,0,0,0,0,0,0,0]))
print(c0.get_syndrome([0,1,1,0,0,0,0,0,0,0]))
print(c0.get_syndrome([1,0,1,0,0,0,0,0,0,0]))
#c2 = STEANE()
c = ConcCode([c0,c1])
dec_sim(c,PROB=[0.13],MONTE=1000,LOG_OUTPUT_SPAN=1);exit()
np.set_printoptions(threshold=10000)
#E = [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1]
E = [0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0 ,0 ,0, 0]
E0 = [0,1,0,0,0,0,0,0,0,0]
p = 0.01
c.set_error_probability(np.array([1-p,p/3,p/3,p/3]),iid=True)
c1.set_error_probability(np.array([1-p,p/3,p/3,p/3]),iid=True)
c0.set_error_probability(np.array([1-p,p/3,p/3,p/3]),iid=True)
result=c.decode(c.get_syndrome(E))
print("T  :",result["T"])
print("L  :",result["L"])
print("LT :",result["LT"])
print("ELT:",E^result["LT"])
print("S  :",c.get_syndrome(E^result["LT"]))
print(c.in_S(E^result["LT"]))
#result0=c0.decode(c0.get_syndrome(E0))
#print(1,c0,result0["LT"])
