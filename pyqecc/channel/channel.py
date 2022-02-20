import numpy as np
from .abstruct import Channel
from ..util import pishifts
class DepolarizingChannel(Channel):
    def __init__(self,p,seed=None):
        super().__init__(seed)
        setp = []
        if type(p) !=list and type(p)!=type(np.array):
            p=[p]
        for pp in p:
            setp.append([1 - pp,pp/3,pp/3,pp/3])
        self._channel_parameter["p"] = np.array(setp)
        self.generate_param()
        
    def channel(self,n,ind=0):
        if n>0 and n!=0:
            self.set_n(n)
        r = np.random.random(n)
        p = 1 - self._channel_parameter["p"][ind][0]
        x_pos = np.where(r <= p / 3)[0]
        z_pos = np.intersect1d(np.where(r < p)[0], np.where(r > 2 * p / 3)[0])
        y_pos = np.intersect1d(np.where(r < 2 * p / 3)[0], np.where(r > p / 3)[0])
        self._channel_output["E"][x_pos] = 1  # X
        self._channel_output["E"][n + z_pos] = 1  # Z
        self._channel_output["E"][y_pos] = 1  # Y
        self._channel_output["E"][n + y_pos] = 1  # Y
        return self._channel_output

class BitFlipChannel(Channel):
    def __init__(self,t,seed=None):
        super().__init__(seed)
        if type(t) !=list and type(t)!=type(np.array):
            t=[t]
        self._channel_parameter["t"] = t
        self.generate_param()

    def channel(self,n,ind=0):
        if n>0 and n!=0:
            self.set_n(n)
        self._channel_output["E"][:self._channel_parameter["t"][ind]] = 1
        np.random.shuffle(self._channel_output["E"])
        return self._channel_output

class PauliChannel(Channel):
    def __init__(self,px,pz=None,seed=None):
        super().__init__(seed)
        if type(px) !=list and type(px)!=type(np.array):
            px=[px]
        self._channel_parameter["px"] = px
        if type(pz) !=list and type(pz)!=type(np.array):
            pz=[pz]
        self._channel_parameter["pz"] = pz
        self.generate_param()

    def channel(self,n,ind=0):
        if n>0 and n!=0:
            self.set_n(n)
        x_pos = np.where(np.random.random(n) <= self._channel_parameter["px"][ind])[0]
        z_pos = np.where(np.random.random(n) <= self._channel_parameter["pz"][ind])[0]
        self._channel_output["E"][x_pos] = 1  # X
        self._channel_output["E"][n + z_pos] = 1  # Z
        return self._channel_output

class GaussianQuantumChannel(Channel):
    def __init__(self,sigma,seed=None,bit_flip=True,phase_flip=True):
        super().__init__(seed)
        self._channel_parameter["sigma"] = sigma
        self._bit_flip = bit_flip
        self._phase_flip = phase_flip
        self.generate_param()

    def set_n(self,n):
        super().set_n(n)
        self._channel_output["DELTA"] = np.zeros(2 * n)

    def channel(self,n=0,ind=0):
        '''
        Return the ""Analog information""
        '''
        if n!=self.n:
            self.set_n(n)
        if self.bit_flip:
            self._channel_output["DELTA"][:self.n] = np.random.normal(scale = self.channel_parameter["sigma"][ind], size  = 2*self.n)[:self.n]
        if self.phase_flip:
            self._channel_output["DELTA"][self.n:] = np.random.normal(scale = self.channel_parameter["sigma"][ind], size  = 2*self.n)[self.n:]
        self._channel_output["E"] *= 0
        
        # 2√π>|E|>√π => error
        delta = pishifts(self.channel_output["DELTA"])
        e_pos = np.where(np.abs(delta)>=np.sqrt(np.pi)/2)[0]
        self._channel_output["E"][e_pos]=1

        return self.channel_output
    

    @property
    def bit_flip(self):
        return self._bit_flip
    @property
    def phase_flip(self):
        return self._phase_flip