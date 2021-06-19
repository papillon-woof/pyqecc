import numpy as np

def symplex_binary_inner_product(a,b):
    n = a.T.shape[0]//2
    z = np.zeros((n,n),dtype='i1')
    i = np.identity(n,dtype='i1')
    Lam = np.c_[np.r_[z,i],np.r_[i,z]]
    return np.mod(np.dot(np.dot(a,Lam),b),2)
