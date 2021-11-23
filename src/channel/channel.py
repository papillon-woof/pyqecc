import numpy as np

def depolarizing_noise(n,p):
    r = np.random.random(n)
    E = np.zeros(2*n,dtype='i1')
    x_pos = np.where(r<=p/3)[0]
    z_pos = np.intersect1d(np.where(r<p)[0], np.where(r>2*p/3)[0])
    y_pos = np.intersect1d(np.where(r<2*p/3)[0], np.where(r>p/3)[0])
    E[x_pos]=1 #X
    E[n+z_pos]=1 #Z
    E[x_pos] = 1;E[n+z_pos] = 1 #Y
    return E

def one_bit_noise(n,p):
    E = np.zeros(2*n,dtype='i1')
    E[np.random.randint(0,n,1)]=1
    return E

def amplitude_damping(n,p):
    pass

