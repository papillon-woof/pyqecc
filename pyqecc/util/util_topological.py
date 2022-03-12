import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

"""
monte = 10000
V=np.zeros((monte,3))
for i in range(monte):
    a = np.random.rand()
    b = np.random.rand()
    V[i][0] = np.sin(a*2*np.pi)*np.cos(b*2*np.pi)
    V[i][1] = np.sin(a*2*np.pi)*np.sin(b*2*np.pi)
    V[i][2] = np.cos(a*2*np.pi)
"""


def plt_V(V):
    X = V[:, 0]
    Y = V[:, 1]
    Z = V[:, 2]
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(X, Y, Z)
    plt.show()


def gen_torus(d1, d2):
    t = np.arange(d1) / d1
    p = np.arange(d2) / d2
    R = 2
    r = 1
    V = []
    E = []
    for i in range(d1):
        for j in range(d2):
            x = R * np.cos(2 * np.pi * t[i]) + r * np.cos(2 * np.pi * p[j]) * np.cos(
                2 * np.pi * t[i]
            )
            y = R * np.sin(2 * np.pi * t[i]) + r * np.cos(2 * np.pi * p[j]) * np.sin(
                2 * np.pi * t[i]
            )
            z = r * np.sin(2 * np.pi * p[j])
            V.append([x, y, z])
    for i in range(d1):
        for j in range(d2):
            if i * d2 + j < ((i + 1) % d1) * d2 + j:
                E.append([i * d2 + j, ((i + 1) % d1) * d2 + j])
            if i * d2 + j < ((i - 1) % d1) * d2 + j:
                E.append([i * d2 + j, ((i - 1) % d1) * d2 + j])
            if i * d2 + j < i * d2 + (j + 1) % d2:
                E.append([i * d2 + j, i * d2 + (j + 1) % d2])
            if i * d2 + j < i * d2 + (j - 1) % d2:
                E.append([i * d2 + j, i * d2 + (j - 1) % d2])
    F = []
    ##面を求める．点の配置順から決定
    for i in range(d1):
        for j in range(d2):
            f = []
            a = i * d2 + j
            b = i * d2 + (j + 1) % d2
            if a > b:
                a, b = b, a
            f.append(E.index([a, b]))
            a = ((i + 1) % d1) * d2 + j
            b = ((i + 1) % d1) * d2 + (j + 1) % d2
            if a > b:
                a, b = b, a
            f.append(E.index([a, b]))
            a = i * d2 + j
            b = ((i + 1) % d1) * d2 + j
            if a > b:
                a, b = b, a
            f.append(E.index([a, b]))
            a = i * d2 + (j + 1) % d2
            b = ((i + 1) % d1) * d2 + (j + 1) % d2
            if a > b:
                a, b = b, a
            f.append(E.index([a, b]))
            F.append(f)

    return np.array(V), np.array(E), np.array(F)


# V,E,F = gen_torus(50,50)
# print(V,E,F)
# print(plt_V(V))
"""
def gen_torus(d):
    theta = (3-d)*np.pi
    origin = np.zeros(3)
    for i in range(d):
        origin -= np.array([np.cos(2*pi*theta),np.cos(2*pi*theta),0])
"""
