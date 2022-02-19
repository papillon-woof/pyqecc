import numpy as np
from .abstruct import Channel
from ..util import dmods
class DepolarizingChannel(Channel):
    def __init__(self,p,seed=None):
        super().__init__(seed)
        setp = []
        for pp in p:
            setp.append([1 - pp,pp/3,pp/3,pp/3])
        self._channel_parameter["p"] = np.array(setp)
        self.generate_param()
        print(self._channel_parameter)

    def channel(self,n,ind=0):
        if n>0 and n!=0:
            self.set_n(n)
        r = np.random.random(n)
        x_pos = np.where(r <= self.p / 3)[0]
        z_pos = np.intersect1d(np.where(r < self.p[ind])[0], np.where(r > 2 * self.p[ind] / 3)[0])
        y_pos = np.intersect1d(np.where(r < 2 * self.p[ind] / 3)[0], np.where(r > self.p[ind] / 3)[0])
        self._channel_output["E"][x_pos] = 1  # X
        self._channel_output["E"][n + z_pos] = 1  # Z
        self._channel_output["E"][y_pos] = 1  # Y
        self._channel_output["E"][n + y_pos] = 1  # Y
        return self._channel_output

class BitFlipChannel(Channel):
    def __init__(self,t,seed=None):
        super().__init__(seed)
        self._channel_parameter["t"] = t
        self.generate_param()

    def channel(self,n,ind=0):
        if n>0 and n!=0:
            self.set_n(n)
        self._channel_output["E"][:self.t[ind]] = 1
        np.random.shuffle(self._channel_output["E"])
        return self._channel_output

class PauliChannel(Channel):
    def __init__(self,px,pz=None,seed=None):
        super().__init__(seed)
        self._channel_parameter["px"] = px
        self._channel_parameter["pz"] = pz
        self.px = px
        self.pz = pz
        self.generate_param()

    def channel(self,n,ind=0):
        if n>0 and n!=0:
            self.set_n(n)
        x_pos = np.where(np.random.random(n) <= self.px[ind])[0]
        z_pos = np.where(np.random.random(n) <= self.pz[ind])[0]
        self._channel_output["E"][x_pos] = 1  # X
        self._channel_output["E"][n + z_pos] = 1  # Z
        return self._channel_output

class GaussianQuantumChannel(Channel):
    def __init__(self,sigma,sigma_z=None,seed=None):
        super().__init__(seed)
        self._channel_parameter["sigma"] = sigma
        self._channel_parameter["sigma_z"] = sigma
        self.sigma_x = sigma
        if sigma_z is None:
            self.sigma_z = sigma
        else:
            self.sigma_z = sigma_z
        self.generate_param()

    def set_n(self,n):
        super().__init__(n)
        self._channel_output["DELTA_M"] = np.zeros(2 * n)

    def channel(self,n=0,ind=0):
        '''
        Return the ""Analog information""
        '''
        if n>0 and n!=0:
            self.set_n(n)
        self._channel_output["DELTA_M"][:self.n] = np.random.rand(self.n)*self.sigma_x[ind]
        self._channel_output["DELTA_M"][self.n:] = np.random.rand(self.n)*self.sigma_z[ind]
        self._channel_output["E"] *= 0
        
        # 2√π>|E|>√π => error 
        delta = dmods(np.abs(self.channel_output["DELTA_M"]),2 * np.sqrt(np.pi))
        e_pos = np.where(delta>np.sqrt(np.pi))[0]
        self._channel_output["E"][e_pos]=1

        return self.channel_output