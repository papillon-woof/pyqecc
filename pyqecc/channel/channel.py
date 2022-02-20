import numpy as np
from .abstruct import Channel
from ..util import pishifts


class DepolarizingChannel(Channel):
    def __init__(self, n, p, seed=None):
        super().__init__(n, seed)
        if not isinstance(p, (list, np.ndarray)):
            p = [p]
        self.param_len = len(p)
        self._channel_parameter["p"] = p

    def channel(self, n=-1, ind=0):
        self._channel_output["E"] *= 0
        if n > 0:
            self.n = n
        r = np.random.random(self.n)
        p = self._channel_parameter["p"][ind]
        x_pos = np.where(r <= p / 3)[0]
        z_pos = np.intersect1d(np.where(r < p)[0], np.where(r > 2 * p / 3)[0])
        y_pos = np.intersect1d(np.where(r < 2 * p / 3)[0], np.where(r > p / 3)[0])
        self._channel_output["E"][x_pos] = 1  # X
        self._channel_output["E"][self.n + z_pos] = 1  # Z
        self._channel_output["E"][y_pos] = 1  # Y
        self._channel_output["E"][self.n + y_pos] = 1  # Y
        return self._channel_output

    def get_param(self, ind):
        p = self.channel_parameter["p"][ind]
        return [1 - p, p / 3, p / 3, p / 3]


class BitFlipChannel(Channel):
    def __init__(self, n, tx, tz, seed=None):
        super().__init__(n, seed)
        if not isinstance(tx, (list, np.ndarray)):
            tx = np.array([tx])
        if not isinstance(tz, (list, np.ndarray)):
            tz = np.array([tz])
        if not isinstance(tx, np.ndarray):
            tx = np.array(tx)
        if not isinstance(tz, np.ndarray):
            tz = np.array(tz)
        self.param_len = min(len(tx), len(tz))
        self._channel_parameter["tx"] = tx
        self._channel_parameter["tz"] = tz
        self.px = self._channel_parameter["tx"] / n
        self.pz = self._channel_parameter["tz"] / n

    def channel(self, n=-1, ind=0):
        self._channel_output["E"] *= 0
        if n > 0:
            self.n = n
            self.px = self._channel_parameter["tx"] / self.n
            self.pz = self._channel_parameter["tz"] / self.n
        (self._channel_output["E"][: self.n])[: self._channel_parameter["tx"][ind]] = 1
        (self._channel_output["E"][self.n :])[: self._channel_parameter["tz"][ind]] = 1
        np.random.shuffle(self._channel_output["E"][: self.n])
        np.random.shuffle(self._channel_output["E"][self.n :])
        return self._channel_output

    def get_param(self, ind):
        return [
            (1 - self.px[ind]) * (1 - self.pz[ind]),
            (1 - self.px[ind]) * self.pz[ind],
            (1 - self.pz[ind]) * self.px[ind],
            self.px[ind] * self.pz[ind],
        ]


class PauliChannel(Channel):
    def __init__(self, n, px, pz, seed=None):
        super().__init__(n, seed)
        if not isinstance(px, (list, np.ndarray)):
            px = [px]
        self._channel_parameter["px"] = px
        if not isinstance(pz, (list, np.ndarray)):
            pz = [pz]
        if not isinstance(px, np.ndarray):
            px = np.array(px)
        if not isinstance(pz, np.ndarray):
            pz = np.array(pz)
        self.param_len = min(len(px), len(pz))
        self._channel_parameter["pz"] = pz

    def channel(self, n=-1, ind=0):
        self._channel_output["E"] *= 0
        if n > 0:
            self.n = n
        x_pos = np.where(
            np.random.random(self.n) <= self._channel_parameter["px"][ind]
        )[0]
        z_pos = np.where(
            np.random.random(self.n) <= self._channel_parameter["pz"][ind]
        )[0]
        self._channel_output["E"][x_pos] = 1  # X
        self._channel_output["E"][self.n + z_pos] = 1  # Z
        return self._channel_output

    def get_param(self, ind):
        px = self._channel_parameter["px"][ind]
        pz = self._channel_parameter["pz"][ind]
        return [
            (1 - px) * (1 - pz),
            (1 - pz) * px,
            (1 - px[ind]) * pz,
            px[ind] * pz[ind],
        ]


class GaussianQuantumChannel(Channel):
    def __init__(self, n, sigma, seed=None, bit_flip=True, phase_flip=True):
        super().__init__(n, seed)
        if not isinstance(sigma, (list, np.ndarray)):
            sigma = np.array([sigma])
        if not isinstance(sigma, (np.ndarray)):
            sigma = np.array(sigma)
        self._channel_parameter["sigma"] = sigma
        self._bit_flip = bit_flip
        self._phase_flip = phase_flip
        self.param_len = len(sigma)

    def channel(self, n=-1, ind=0):
        """
        Return the ""Analog information""
        """
        if n > 0:
            self.n = n
        if self.bit_flip:
            self._channel_output["DELTA"][: self.n] = np.random.normal(
                scale=self.channel_parameter["sigma"][ind], size=2 * self.n
            )[: self.n]
        if self.phase_flip:
            self._channel_output["DELTA"][self.n :] = np.random.normal(
                scale=self.channel_parameter["sigma"][ind], size=2 * self.n
            )[self.n :]
        self._channel_output["E"] *= 0
        # 2√π>|E|>√π => error
        delta = pishifts(self.channel_output["DELTA"])
        e_pos = np.where(np.abs(delta) >= np.sqrt(np.pi) / 2)[0]
        self._channel_output["E"][e_pos] = 1

        return self.channel_output

    def get_param(self, ind):
        return self._channel_parameter["sigma"][ind]

    @property
    def bit_flip(self):
        return self._bit_flip

    @property
    def phase_flip(self):
        return self._phase_flip

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, n):
        if n <= 0:
            raise ValueError("PHYSICAL_QUBIT is more than 0")
        self._n = n
        self._channel_output["E"] = np.zeros(2 * n, dtype="i1")
        self._channel_output["DELTA"] = np.zeros(2 * n)
