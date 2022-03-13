import numpy as np
import networkx as nw
from .abstruct import *
from .stabilizer import *
from ..util import *
import csv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import datetime
from datetime import datetime

def plt_V(V):
    if type(V) == list:
        V = np.array(V)
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
    edges = []
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
                edges.append([i * d2 + j, ((i + 1) % d1) * d2 + j])
            if i * d2 + j < ((i - 1) % d1) * d2 + j:
                edges.append([i * d2 + j, ((i - 1) % d1) * d2 + j])
            if i * d2 + j < i * d2 + (j + 1) % d2:
                edges.append([i * d2 + j, i * d2 + (j + 1) % d2])
            if i * d2 + j < i * d2 + (j - 1) % d2:
                edges.append([i * d2 + j, i * d2 + (j - 1) % d2])
    F = []
    ##面を求める．点の配置順から決定
    for i in range(d1):
        for j in range(d2):
            f = []
            a = i * d2 + j
            b = i * d2 + (j + 1) % d2
            if a > b:
                a, b = b, a
            f.append(edges.index([a, b]))
            a = ((i + 1) % d1) * d2 + j
            b = ((i + 1) % d1) * d2 + (j + 1) % d2
            if a > b:
                a, b = b, a
            f.append(edges.index([a, b]))
            a = i * d2 + j
            b = ((i + 1) % d1) * d2 + j
            if a > b:
                a, b = b, a
            f.append(edges.index([a, b]))
            a = i * d2 + (j + 1) % d2
            b = ((i + 1) % d1) * d2 + (j + 1) % d2
            if a > b:
                a, b = b, a
            f.append(edges.index([a, b]))
            F.append(f)

    # edgesにインデックスをつける処理
    return (
        np.r_[[np.arange(len(V))], np.array(V).T].T,
        np.r_[[np.arange(len(edges))], np.array(edges).T].T,
        np.array(F),
    )

class Surface_old(SC):
    # 点，辺，面の作りは，以下で規定されている．
    # 点は，座標および添え字，辺は点の組み合わせ，面は辺の番号の添え字
    # 双対点は重心および，面の番号．双対辺は添え字は対応している．
    _name = "surface code"

    def __init__(self, fname, H="obj", T=None, L=None, P=None, iid=True):
        # H is replesentation for graph into adjacent matrix:
        self._fname = fname
        if H == "obj":
            self.read_obj()
        elif H == "csv":
            self.read_txt()
            self._H = H
        elif H == "H":
            self._H = H

        self.enc_circuit = None
        self.dec_circuit = None
        self._P = self.set_P(P)
        self.ML_decoding_qubit_limit = 15

    def read_obj(self):
        pass

    def read_txt(self):
        # 形式
        # 頂点: 点の座標，行番号が点番号: float[3]
        # 辺: 点の集合 edge 行番号が辺番号: int[2]
        # 面: 辺の集合 face 行番号が面番号: int[d]
        # 現状3次元まで
        # 12 15 17
        # 1 1.0 2.1 3,1
        # 2 1.1 2.3 4.5
        # ...
        # 1 2 5
        # 1 3
        # 1 2
        # 1 2 4 5..
        # 1 3 5 6..
        # ...

        with open(self.fname, newline="") as f:
            reader = [row for row in csv.reader(f)]
            self._num_v, self._num_e, self._num_f = [int(i) for i in reader[0]]
            self._V = [[] for i in range(self.num_v)]
            self._edges = [[] for i in range(self.num_e)]
            self._F = [[] for i in range(self.num_f)]
            self._VV = [[] for i in range(self.num_f)]
            self._EE = [[] for i in range(self.num_e)]
            self._FF = [[] for i in range(self.num_v)]

            for i in range(1, self.num_v + 1):
                v = [float(data) for data in reader[i]]
                self._V[int(v[0])] = v[1:]
            for i in range(self.num_e):
                e = [int(data) for data in reader[i + self.num_v + 1]]
                self._E[int(e[0])] = e[1:]
            for i in range(self.num_f):
                f = [int(data) for data in reader[1 + i + self.num_e + self.num_v]]
                self._F[i] = f
        print(1, self.V)
        print(2, self.E)
        print(3, self.F)

        # 双対辺
        # self._EE = [[] for i in range(self.num_e)]
        for i in range(self.num_f):
            for j in range(len(self.F[i])):
                self._EE[self.F[i][j]].append(i)
        print(self._EE)

        # 双対点(面の重心を点とする)
        for i in range(self.num_f):
            v_all = []
            e_all = []
            for e in self.F[i]:
                e_all.append(self.edges[e][0])
                e_all.append(self.edges[e][1])
            e_all = list(set(e_all))
            for e in e_all:
                v_all.append(self.V[e])
            v_xyz = np.zeros(3)
            for v in v_all:
                v_xyz[0] += v[0]
                v_xyz[1] += v[1]
                v_xyz[2] += v[2]
            v_xyz /= len(v_all)
            self._VV[i] = v_xyz.tolist()
        print(self._VV)
        # plt_V(self._VV)

        # 双対面
        for i in range(self.num_e):
            self._FF[self.edges[i][0]].append(i)
            self._FF[self.edges[i][1]].append(i)
        print(self._FF)

        # パリティ検査行列作成
        self._n = self.num_e
        self._k = self.num_e - ((self.num_v - 1) + (self.num_f - 1))
        self._R = self.k / self.n
        self._H = np.zeros((self.n - self.k, 2 * self.n), dtype="i1")
        for i in range(self.num_f - 1):
            print(self.F[i], self.num_f - 1)
            self._H[i][self.F[i]] = 1
        for i in range(self.num_v - 1):
            self._H[self.num_f - 1 + i][np.array(self.FF[i]) + self.num_e] = 1
        print(self.H)

    def decode(self, s):
        return self.mwpm(s)

    def matching(self, syndromes):
        graph_instance = nw.Graph()
        edges = self.get_qubit_distances(syndromes, self.code.size)
        for v0, v1, weight in edges:
            graph_instance.add_edge(v0, v1, weight=-weight)
        return nw.algorithms.matching.max_weight_matching(graph_instance, maxcardinality=10)

    def correct_matching(self, syndromes, matching, **kwargs):
        weight = 0
        for i, j in matching:
            weight += self._correct_matched_qubits(syndromes[i], syndromes[j])
        return weight

    def get_qubit_distances(self,syndrome):
        '''
        calculated between qubits. wo measurement err
        s: v
        重みを計算する必要がある．
        toricでは変更
        '''
        edges = []
        for i, s in enumerate(syndrome):
            if s == 1:
                x, y = self._FE[i]
            for j, ss in enumerate(syndrome[i + 1:]):
                xx, yy = self._FE[j]
                wx = np.abs(x - xx)
                wy = np.abs(y - yy)
                weight = wx + wy 
                edges.append([i, j + i + 1, weight])
        print(edges)
        return edges

    def mwpm(self, s):
        x = self.correct_matching(s[:len(self.nk)//2], self.matching(s[:len(self.nk)//2]))
        z = self.correct_matching(s[len(self.nk)//2:], self.matching(s[len(self.nk)//2:]))
        return x+z

    def set_P(self, P):
        self._P = P

    @property
    def num_v(self):
        return self._num_v

    @property
    def num_e(self):
        return self._num_e

    @property
    def num_f(self):
        return self._num_f

    @property
    def V(self):
        return self._V

    @property
    def E(self):
        return self._E

    @property
    def F(self):
        return self._F

    @property
    def VV(self):
        return self._VV

    @property
    def EE(self):
        return self._EE

    @property
    def FF(self):
        return self._FF

    @property
    def fname(self):
        return self._fname

class Surface(SC):
    # 点，辺，面の作りは，以下で規定されている．
    # 点は，座標および添え字，辺は点の組み合わせ，面は辺の番号の添え字
    # 双対点は重心および，面の番号．双対辺は添え字は対応している．
    _name = "surface code"

    def __init__(self, fname, H="obj", T=None, L=None, P=None, iid=True):
        # H is replesentation for graph into adjacent matrix:
        self._fname = fname
        if H == "obj":
            self.read_obj()
        elif H == "csv":
            self.read_txt()
            self._H = H
        elif H == "H":
            self._H = H

        self.enc_circuit = None
        self.dec_circuit = None
        self.ML_decoding_qubit_limit = 15

    def read_txt(self):
        # 形式
        # 頂点: 点の座標，行番号が点番号: float[3]
        # 辺: 点の集合 edge 行番号が辺番号: int[2]
        # 面: 辺の集合 face 行番号が面番号: int[d]
        # 現状3次元まで
        # 12 15 17
        # 1 1.0 2.1 3,1
        # 2 1.1 2.3 4.5
        # ...
        # 1 2 5
        # 1 3
        # 1 2
        # 1 2 4 5..
        # 1 3 5 6..
        # ...
        with open(self.fname, newline="") as f:
            reader = [row for row in csv.reader(f)]
            self._num_v, self._num_e, self._num_f = [int(i) for i in reader[0]]
            self._V = [[],[[] for i in range(self.num_f)]] #双対点は面の個数と同等
            self._E = [[],[[] for i in range(self.num_e)]]
            self._F = [[],[]]
            count = 0
            print(len(reader))
            for i in range(1, self.num_v + 1):
                v = [float(data) for data in reader[i]]
                self._V[0].append(v)
                count+=1
            for i in range(self.num_e):
                e = [int(data) for data in reader[count]]
                self._E[0].append(e)
                count+=1
            for i in range(self.num_f):
                f = [int(data) for data in reader[count]]
                self._F[0].append(f)
                count+=1
            for i in range(1, self.num_v + 1):
                v = [float(data) for data in reader[count]]
                self._V[0].append(v)
                count+=1
            print(count)
            for i in range(self.num_e):
                e = [int(data) for data in reader[count]]
                self._E[0].append(e)
                count+=1
            for i in range(self.num_f):
                f = [int(data) for data in reader[count]]
                self._F[0].append(f)
                count+=1
        
        # Generate a parity chack matrix 
        self._n = self.num_e
        self._k = self.num_e - 2*self.num_f
        self._nk = self.n - self.k
        self._R = self.k / self.n
        self._H = np.zeros((self.nk, 2 * self.n), dtype="i1")
        for i in range(self.num_f):
            self._H[i][self.dual_faces[i]] = 1
        for i in range(self.num_f):
            self._H[self.nk//2 + i][self.n + np.array(self.dual_faces[i])] = 1

    def decode(self, s):
        return self.mwpm(s)

    def matching(self, syndromes,val):
        graph_instance = nw.Graph()
        edges = self.get_qubit_distances(syndromes, val)
        for v0, v1, weight in edges:
            graph_instance.add_edge(v0, v1, weight=-weight)
        return nw.algorithms.matching.max_weight_matching(graph_instance, maxcardinality=10)

    def correct_matching(self, syndromes, matching, **kwargs):
        weight = 0
        for i, j in matching:
            weight += self._correct_matched_qubits(syndromes[i], syndromes[j])
        return weight

    def get_qubit_distances(self,syndrome,val):
        '''
        calculated between qubits. wo measurement err
        s: v
        重みを計算する必要がある．
        toricでは変更
        '''
        edges = []
        vertexs = self.vertexs if val == 0 else self.dual_vertexs
        for i, f in enumerate(syndrome):
            if -1 not in vertexs[i]:
                x,y,z = vertexs[i]
            if s == 1:
                for j, ss in enumerate(syndrome[i + 1:]):
                    if ss == 1:
                        xx,yy,zz = self._vertex[(val+1)%2][j]
                        weight = np.abs(xx - x) + np.abs(yy - y) + np.abs(zz - z)
                        edges.append([i, j + i + 1, weight])
        return edges

    def mwpm(self, s):
        x = self.correct_matching(s[:self.nk//2], self.matching(s[:self.nk//2],0))
        z = self.correct_matching(s[self.nk//2:], self.matching(s[self.nk//2:],1))
        return x+z

    def set_P(self, P):
        self._P = P

    @property
    def num_v(self):
        return self._num_v

    @property
    def num_e(self):
        return self._num_e

    @property
    def num_f(self):
        return self._num_f

    @property
    def vertexs(self):
        return self._vertexs

    @property
    def edges(self):
        return self._edges

    @property
    def faces(self):
        return self._faces

    @property
    def dual_vertexs(self):
        return self._dual_vertexs

    @property
    def dual_edges(self):
        return self._dual_edges

    @property
    def dual_faces(self):
        return self._dual_faces

    @property
    def fname(self):
        return self._fname


class Toric(Surface):
    def __init__(self, d1, d2, P=None):
        self._d1 = d1
        self._d2 = d2
        self._V = [[] for i in range(2)]
        self._E = [[] for i in range(2)]
        self._F = [[] for i in range(2)]
        self._V[0], self._E[0], self._F[0] = gen_torus(self.d1, self.d2)
        self._num_v = len(self.V[0])
        self._num_e = len(self.E[0])
        self._num_f = len(self.F[0])
        if not os.path.isdir('pyqecc/qecc/topological_data'):
            os.makedirs('pyqecc/qecc/topological_data')
        self._fname = 'pyqecc/qecc/topological_data'+'/toric_'+str(d1)+'_'+str(d2)+'_'+datetime.now().strftime('%Y%m%d%H%M%S')+'.csv'
        with open(self.fname, 'w',newline="") as f:
            writer = csv.writer(f)
            writer.writerow([self.num_v, self.num_e, self.num_f])
            for i in range(self.num_v):
                writer.writerow(self.V[i])
            for i in range(self.num_e):
                writer.writerow(self.E[i])
            for i in range(self.num_f):
                writer.writerow(self.F[i])
        self.read_txt()
        self.set_T()
        self.set_L()
        self.enc_circuit = None
        self.dec_circuit = None
        self._n = self.H.shape[0]
        self._nk = self.H.shape[1]
        self._k = self.n - self.nk

    def set_T(self, T=None):
        pass

    def set_L(self, L=None):
        pass

    def decode(self, syndrome,mode="mwpm"):
        super().decode(syndrome)

    @property
    def d1(self):
        return self._d1

    @property
    def d2(self):
        return self._d2

    @property
    def fname(self):
        return self._fname

    def __str__(self):
        output = ""
        output += "codename        :" + str(self.name) + "\n"
        output += "n               : " + str(self.n) + "\n"
        output += "k               : " + str(self.k) + "\n"
        output += "R               : " + str(self.R) + "\n"
        return output


def gen_plane(d1,d2):
    vertexs = []
    edges = []
    faces = []
    dual_vertexs = []
    dual_edges = []
    dual_faces = []
    '''
    - - - -
     | | |
    - - - - 
     | | |
    - - - -
     | | |
    - - - -  
    '''
    # Add the vertexes for planner
    for i in range(d1):
        for j in range(d2+1):
            vertexs.append([i,j,0])
    
    # Add the edges for planner    
    for i in range(d1):
        for j in range(d2):
            edges.append([vertexs.index([i,j,0]),vertexs.index([i,j+1,0])])
    for i in range(d1-1):
        for j in range(1,d2):
            edges.append([vertexs.index([i,j,0]),vertexs.index([i+1,j,0])])

    # Add the faces for planner
    for i in range(d1-1):
        for j in range(d2):    
            f = []
            a = vertexs.index([i,j,0])
            b = vertexs.index([i,j+1,0])
            if a>b:
                a,b = b,a
            
            f.append(edges.index([a, b]))
            a = vertexs.index([i+1,j,0])
            b = vertexs.index([i+1,j+1,0])
            if a>b:
                a,b = b,a
            
            f.append(edges.index([a,b]))
            if j != 0:
                a = vertexs.index([i+1,j,0])
                b = vertexs.index([i,j,0])
                if a>b:
                    a,b = b,a
            
                f.append(edges.index([a,b]))
            if j != d2-1:
                a = vertexs.index([i,j+1,0])
                b = vertexs.index([i+1,j+1,0])
                if a>b:
                    a,b = b,a
                f.append(edges.index([a,b]))
            faces.append(f)
    
    # Mask the vertexes for planner on boundary
    for i in range(d1):
        for j in range(d2+1):
            if j == 0 or j == d2:
                vertexs[i*(d2+1) + j] = [-1]

    # Add the vertexes for planner
    for i in range(d1+1):
        for j in range(d2):
            dual_vertexs.append([i,j,0])
    
    
    # Add the edges for planner    
    for i in range(d1):
        for j in range(d2):
            dual_edges.append([dual_vertexs.index([i,j,0]),dual_vertexs.index([i+1,j,0])])
    for i in range(1,d1):
        for j in range(d2-1):
            dual_edges.append([dual_vertexs.index([i,j,0]),dual_vertexs.index([i,j+1,0])])
    print(vertexs, edges, faces, dual_vertexs, dual_edges, dual_faces,)
    # Add the faces for planner
    for i in range(d1):
        for j in range(d2-1):    
            f = []
            if i != 0:
                a = dual_vertexs.index([i,j,0])
                b = dual_vertexs.index([i,j+1,0])
                if a>b:
                    a,b = b,a
                f.append(dual_edges.index([a, b]))
            if i != d1-1:
                a = dual_vertexs.index([i+1,j,0])
                b = dual_vertexs.index([i+1,j+1,0])
                if a>b:
                    a,b = b,a
                f.append(dual_edges.index([a,b]))
            a = dual_vertexs.index([i+1,j,0])
            b = dual_vertexs.index([i,j,0])
            if a>b:
                a,b = b,a
            f.append(dual_edges.index([a,b]))
            
            a = dual_vertexs.index([i,j+1,0])
            b = dual_vertexs.index([i+1,j+1,0])
            if a>b:
                a,b = b,a
            f.append(dual_edges.index([a,b]))
            dual_faces.append(f)
    
    # Mask the vertexes for planner on boundary
    for i in range(d1+1):
        for j in range(d2):
            if i == 0 or i == d1:
                dual_vertexs[i*d2 + j] = [-1]

    return vertexs, edges, faces, dual_vertexs, dual_edges, dual_faces,

class Planner(Surface):
    def __init__(self, d1, d2, P=None):
        self._d1 = d1
        self._d2 = d2
        self._vertexs = []
        self._edges = []
        self._faces = []
        self._dual_vertexs = []
        self._dual_edges = []
        self._dual_faces = []
        self._vertexs, self._edges, self._faces,self._dual_vertexs, self._dual_edges, self._dual_faces = gen_plane(self.d1, self.d2)
        self._num_v = len(self.vertexs)
        self._num_e = len(self.edges)
        self._num_f = len(self.faces)
        if not os.path.isdir('pyqecc/qecc/topological_data'):
            os.makedirs('pyqecc/qecc/topological_data')
        self._fname = 'pyqecc/qecc/topological_data'+'/toric_'+str(d1)+'_'+str(d2)+'_'+datetime.now().strftime('%Y%m%d%H%M%S')+'.csv'
        with open(self.fname, 'w',newline="") as f:
            writer = csv.writer(f)
            writer.writerow([self.num_v, self.num_e, self.num_f])
            print(self.num_e,self.num_v,self.num_f,)
            for i in range(self.num_v):
                writer.writerow(self.vertexs[i])
            for i in range(self.num_e):
                writer.writerow(self.edges[i])
            for i in range(self.num_f):
                writer.writerow(self.faces[i])
            for i in range(self.num_v):
                writer.writerow(self.dual_vertexs[i])
            for i in range(self.num_e):
                writer.writerow(self.dual_edges[i])
            for i in range(self.num_f):
                writer.writerow(self.dual_faces[i])
        self.read_txt()
        self.enc_circuit = None
        self.dec_circuit = None
        self._n = self.H.shape[0]
        self._nk = self.H.shape[1]
        self._k = self.n - self.nk

    def set_T(self, T=None):
        pass

    def set_L(self, L=None):
        pass

    def decode(self, syndrome,mode="mwpm"):
        super().decode(syndrome)

    @property
    def d1(self):
        return self._d1

    @property
    def d2(self):
        return self._d2