import numpy as np

def STEANE():
    H = np.array([[0,0,0,1,1,1,1,0,0,0,0,0,0,0],[0,1,1,0,0,1,1,0,0,0,0,0,0,0],[1,0,1,0,1,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,1,1,1],[0,0,0,0,0,0,0,0,1,1,0,0,1,1],[0,0,0,0,0,0,0,1,0,1,0,1,0,1],dtype='i1')
    return SC(7,4,H)
