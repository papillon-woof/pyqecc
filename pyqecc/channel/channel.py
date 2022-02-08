from abc import ABCMeta, abstractmethod
import numpy as np

def depolarizing(n,p):
    r = np.random.random(n)
    E = np.zeros(2 * n,dtype='i1')
    x_pos = np.where(r<=p/3)[0]
    z_pos = np.intersect1d(np.where(r<p)[0], np.where(r>2*p/3)[0])
    y_pos = np.intersect1d(np.where(r<2*p/3)[0], np.where(r>p/3)[0])
    E[x_pos]=1 #X
    E[n+z_pos]=1 #Z
    E[y_pos] = 1;E[n+y_pos] = 1 #Y
    return E

def one_bit_flip(n,p):
    E = np.zeros(2*n,dtype='i1')
    E[np.random.randint(0,n,1)]=1
    return E

def amplitude_damping(n,p):
    pass

def channel(n,p,CHANNEL_MODEL="DEPOLARIZING"):
    if CHANNEL_MODEL=="DEPOLARIZING":
        return depolarizing(n,p)
    elif CHANNEL_MODEL=="ONE_BIT_FLIP":
        return depolarizing(n,p)
    elif CHANNEL_MODEL=="AMPLITUDE_DUMPING":
        ValueError("Not implimentation")
        #return amplitude_damping(n,p)

def gaussian_quantum_channel(n,p):
    pass

def photon_loss_channel(n,p):
    pass
