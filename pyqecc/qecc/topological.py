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

    def decode(self, s):
        return self.mwpm(s)

    def matching(self, syndromes,val):
        # minimum weight perfect matching
        graph_instance = nw.Graph()
        edges = self.get_qubit_distances(syndromes, val)
        for v0, v1, weight in edges:
            graph_instance.add_edge(v0, v1, weight=-weight)
        return nw.algorithms.matching.max_weight_matching(graph_instance, maxcardinality=10)

    def correct_matching(self, matching, val, **kwargs):
        # get_error
        qubit = set()
        for i, j in matching:
            qubit^=self._correct_matched_qubits(i, j, val)
        ret = np.zeros(2*self.n,dtype='i1')
        ret[list(qubit)] = 1
        return ret

    def get_qubit_distances(self,syndrome,val):
        '''
        calculated between qubits. wo measurement err
        s: v
        '''
        edges = []
        for s_ind, s in enumerate(syndrome):
            if s != 1:
                continue
            f = np.where(self.H[self.nk//2*val+s_ind]==1)[0]
            f.sort()
            f = f.tolist()
            for i in range(self.d[val]-1):
                if f in self.faces[val][i]:
                    x,y = i,self.faces[val][i].index(f)
                    break
            for ss_ind, ss in enumerate(syndrome[s_ind + 1:]):
                if ss != 1:
                    continue
                f = np.where(self.H[self.nk//2*val+s_ind+ss_ind+1]==1)[0]
                f.sort()
                for j in range(self.d[val]-1):
                    if f.tolist() in self.faces[val][j]:
                        xx,yy = j,self.faces[val][j].index(f.tolist())
                        weight = np.abs(xx-x)+np.abs(yy-y)
                        edges.append([s_ind, s_ind + ss_ind + 1,weight])
                        break
        return edges

    def _correct_matched_qubits(self, s0, s1, val):
        """Flips the values of edges between two matched qubits by doing a walk in between."""   
        x0,y0 = s0%(self.d[0]), s0//(self.d[0])
        x1,y1 = s1%(self.d[0]), s1//(self.d[0])
        dx,dy = np.sign(x1 - x0),np.sign(y1 - y0)
        qubit = set()
        while(y0!=y1):
            qubit^=set(self.faces[val][y0+dy][x0])&set(self.faces[val][y0][x0])
            y0+=dy
        while(x0!=x1):
            qubit^=set(self.faces[val][y0][x0+dx])&set(self.faces[val][y0][x0])
            x0+=dx
        return qubit

    def mwpm(self, s):
        x = self.correct_matching(self.matching(s[:self.nk//2],0),0)
        z = self.correct_matching(self.matching(s[self.nk//2:],1),1)
        return x+z

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
    def fname(self):
        return self._fname


class Toric(Surface):
    def __init__(self, d1, d2, P=None):
        self._d1 = d1
        self._d2 = d2
        self._V = [[] for i in range(2)]
        self._E = [[] for i in range(2)]
        self._F = [[] for i in range(2)]
        self._V[0], self._E[0], self._F[0] = gen_torus(self.d[0], self.d[1])
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

    # original d1 x (d2 + 1)
    # dual: (d1 - 1) x d2

    vertexs = []
    edges = []
    faces = [[[] for i in range(d2)] for j in range(d1-1)]
    dual_vertexs = []
    dual_edges = []
    dual_faces = [[[] for i in range(d2 - 1)] for j in range(d1)]

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
            faces[i][j] = f
    
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
            dual_faces[i][j] = f
    
    # Mask the vertexes for planner on boundary
    for i in range(d1+1):
        for j in range(d2):
            if i == 0 or i == d1:
                dual_vertexs[i*d2 + j] = [-1]

    # faces <-> dual_vertexs
    graph_component = {"vertexs":vertexs, "edges":edges, "faces":faces, "dual_vertexs":dual_vertexs, "dual_edges":dual_edges, "dual_faces":dual_faces,}
    faces_coordinate = []
    dual_faces_coordinate = []

    return vertexs, edges, faces, dual_vertexs, dual_edges, dual_faces,

class Planner(Surface):
    def __init__(self, d1, d2, P=None):
        self._d = [d1,d2]
        self._vertexs = [[] for _ in range(2)]
        self._edges = [[] for _ in range(2)]
        self._faces = [[] for _ in range(2)]
        self._vertexs[0], self._edges[0], self._faces[0], self._vertexs[1], self._edges[1], self._faces[1] = gen_plane(self.d[0], self.d[1])
        self._num_v = [len(self.vertexs[0]),len(self.vertexs[1])]
        self._num_e = [len(self.edges[0]),len(self.edges[1])]
        self._num_f = [(self.d[0] - 1) * self.d[1],(self.d[1] - 1) * self.d[0]]
        if not os.path.isdir('pyqecc/qecc/topological_data'):
            os.makedirs('pyqecc/qecc/topological_data')
        self._fname = 'pyqecc/qecc/topological_data'+'/toric_'+str(self.d[0])+'_'+str(self.d[1])+'_'+datetime.now().strftime('%Y%m%d%H%M%S')+'.csv'
        with open(self.fname, 'w',newline="") as f:
            writer = csv.writer(f)
            writer.writerow([self.num_v[0], self.num_e[0], self.num_f[0]])
            writer.writerow([self.num_v[1], self.num_e[1], self.num_f[1]])
            for i in range(2):
                for j in range(self.num_v[i]):
                    writer.writerow(self.vertexs[i][j])
                for j in range(self.num_e[i]):
                    writer.writerow(self.edges[i][j])
                for j in range(self.num_f[i]):
                    writer.writerow(self.faces[i][j//(self.d[i] - 1)][j%(self.d[i] - 1)])
        self.read_txt()
        self.enc_circuit = None
        self.dec_circuit = None
        self._n = self.H.shape[1]//2
        self._nk = self.H.shape[0]
        self._k = self.n - self.nk

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
            vertexs = [[] for _ in range(2)]
            edges = [[] for _ in range(2)]
            faces = [[] for _ in range(2)]
            count = 2
            for i in range(2):
                for j in range(self.num_v[i]):
                    v = [float(data) for data in reader[count]]
                    vertexs[i].append(v)
                    count+=1
                for j in range(self.num_e[i]):
                    e = [int(data) for data in reader[count]]
                    edges[i].append(e)
                    count+=1
                for j in range(self.d[(i+1)%2]):
                    f = []
                    for k in range(self.d[i]-1):
                        f.append([int(data) for data in reader[count]])
                        count+=1
                    faces[i].append(f)
        self._vertexs = vertexs
        self._edges = edges
        self._faces = faces
        print(faces[1])
        # Generate a parity chack matrix 
        self._n = self.num_e[0]
        self._k = self.num_e[0] - sum(self.num_f)
        self._nk = self.n - self.k
        self._R = self.k / self.n
        self._H = np.zeros((self.nk, 2 * self.n), dtype="i1")

        for i in range(2):
            for j in range(self.num_f[i]):
                self._H[i*self.num_f[i]+j][self.faces[i][j//(self.d[i] - 1)][j%(self.d[i] - 1)]] = 1
    
    def decode(self, syndrome,mode="mwpm"):
        return super().decode(syndrome)

    @property
    def d(self):
        return self._d