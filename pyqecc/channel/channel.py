from .abstruct import *
import numpy as np

class DepolarizingChannel(Channel):
    def __init__(self,p,seed=None):
        self.p = p
        super().__init__(seed)

    def channel(self,n):
        r = np.random.random(n)
        E = np.zeros(2 * n, dtype="i1")
        x_pos = np.where(r <= self.p / 3)[0]
        z_pos = np.intersect1d(np.where(r < self.p)[0], np.where(r > 2 * self.p / 3)[0])
        y_pos = np.intersect1d(np.where(r < 2 * self.p / 3)[0], np.where(r > self.p / 3)[0])
        E[x_pos] = 1  # X
        E[n + z_pos] = 1  # Z
        E[y_pos] = 1
        E[n + y_pos] = 1  # Y
        return E

class TBitFlipChannel(Channel):
    def __init__(self,t,seed=None):
        self.t = t
        super().__init__(seed)

    def channel(self,n):
        E = np.zeros(2 * n, dtype="i1")
        E[:self.t] = 1
        np.random.shuffle(E)
        return E

class PauliChannel(Channel):
    def __init__(self,px,pz,seed=None):
        self.px = px
        self.pz = pz
        super().__init__(seed)

    def channel(self,n):
        E = np.zeros(2 * n, dtype="i1")
        x_pos = np.where(np.random.random(n) <= self.px)[0]
        z_pos = np.where(np.random.random(n) <= self.pz)[0]
        E[x_pos] = 1  # X
        E[n + z_pos] = 1  # Z
        return E