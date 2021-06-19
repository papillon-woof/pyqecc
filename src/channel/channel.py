import numpy as np

def depolarizing_noise(n,p):
    e = np.zeros(2*n,dtype='i1')
    return e

print(depolarizing_noise(10,0.1))
