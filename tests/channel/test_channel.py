import numpy as np
from pyqecc.channel import DepolarizingChannel, BitFlipChannel, PauliChannel


def test_depolarizing():
    n = 100000
    for p in [0.01, 0.1, 0.2]:
        myChannel = DepolarizingChannel(n, p=p, seed=0)
        e = myChannel.channel(n)["E"]
        cx = 0
        cy = 0
        cz = 0
        for i in range(n):
            if e[i] == 1 and e[i + n] == 1:
                cy += 1
            elif e[i] == 1:
                cx += 1
            elif e[i + n] == 1:
                cz += 1
        assert np.abs(cy / n - p / 3) < 0.003
        assert np.abs(cx / n - p / 3) < 0.003
        assert np.abs(cz / n - p / 3) < 0.003


def test_t_bitflip():
    n = 1000
    for t in range(10):
        myChannel = BitFlipChannel(n, tx=t, tz=t, seed=0)
        for i in range(30):
            e = myChannel.channel(n)["E"]
            assert sum(e) == 2 * t


def test_pauli():
    n = 100000
    for p in [0.01, 0.1, 0.2, 0.3, 0.02, 0.4]:
        myChannel = PauliChannel(n, px=p, pz=p / 2, seed=0)
        e = myChannel.channel(ind=0)["E"]
        cx = 0
        cz = 0
        for i in range(n):
            if e[i] == 1:
                cx += 1
            if e[i + n] == 1:
                cz += 1
        assert np.abs(cx / n - p) < 0.003
        assert np.abs(cz / n - p / 2) < 0.003
