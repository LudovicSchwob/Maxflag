#codé à FPSAC 2025
#cf. conversation avec Germain Pouillot

import matplotlib.pyplot as plt

from collections import defaultdict
def RhombusTilingsGraph(l):
    """
    l = [[1,2,3],[1,4,5],[4,2,6],[6,3,5]]
    """
    D = defaultdict(list)
    O = []
    for i,p in enumerate(l):
        for k in p:
            D[k].append(i)
        O.append({k:j for j,k in enumerate(p)})
    l = tuple(tuple(t) for t in l)
    L = [l]
    LG = {}
    C3 = graphs.CycleGraph(3)
    while len(L)>0:
        L2 = set()
        for l in L:
            lv = []
            m = []
            for i,p in enumerate(l):
                for j in range(len(p)-1):
                    if O[i][p[j]] < O[i][p[j+1]]:
                        m.append((p[j],p[j+1]))
            G = Graph(m)
            for g in G.subgraph_search_iterator(C3, return_graphs=False):
                if g[0]<g[1]<g[2]:
                    l2 = [list(t) for t in l]
                    for k in g:
                        for i in D[k]:
                            j = l[i].index(k)
                            if j+1<len(l[i]) and l[i][j+1] in g:
                                l2[i][j],l2[i][j+1] = l2[i][j+1],l2[i][j]
                    l2 = tuple(tuple(t) for t in l2)
                    lv.append(l2)
                    L2.add(l2)
            LG[l] = lv
        L = L2
    return LG

def RhombusTiling(n):
    l,k = [],1
    for i in range(n):
        m = list(range(k,k+n-i-1))
        for j in range(i):
            m.append(l[-j-1][-i])
        l.append(m)
        k += n-i-1
    return l

def ReducedWord_to_RhombusTiling(n,w):
    """
    w = mot réduit pour w0
    par ex. w = [1,2,1,3,2,1,4,3,2,1]
    """
    if len(w)!=n*(n-1)//2:
        raise Exception("w n'a pas la bonne longueur")
    p,labels = list(range(1,n+1)),[]
    for j in w:
        labels.append((p[j]-1)*(p[j]-2)//2+p[j-1])
        p[j-1],p[j] = p[j],p[j-1]
    T = []
    for k in range(n):
        state = k
        t = []
        for i,j in enumerate(w):
            if j == state:
                state -= 1
                t.append(labels[i])
            elif j == state + 1:
                state += 1
                t.append(labels[i])
        T.append(t)
    return T

#needs 'Middle orders.sage'

def DistributiveSortingOrders(n):
    W = WeylGroup(['A',n-1],prefix = 's')
    LL = WSortingOrders(W)
    element_colors = {'green':[],'red':[]}
    for L in LL:
        T = tuple(tuple(t) for t in ReducedWord_to_RhombusTiling(n,LL[L]))
        if L.is_distributive():
            element_colors['green'].append(T)
        else:
            element_colors['red'].append(T)
    minw = []
    for k in range(1,n):
        minw.extend(list(range(k,0,-1)))
    return Poset(RhombusTilingsGraph(ReducedWord_to_RhombusTiling(n,minw))), element_colors

def RandomReducedWord(n):
    p,w = list(range(1,n+1)),[]
    for _ in range(n*(n-1)//2):
        asc = [i for i in range(1,n) if p[i-1]<p[i]]
        r = asc[randint(0,len(asc)-1)]
        w.append(r)
        p[r],p[r-1] = p[r-1],p[r]
    return w

class W0ReducedWord():
    def __init__(self, n, w):
        if len(w) != n*(n-1)//2:
            raise Exception('w must have length n*(n-1)/2')
        self.N = n
        self.main = w
    def __str__(self):
        return str(self.main)
    def __repr__(self):
        return str(self.main)
    def __eq__(self, other):
        if type(other)!=W0ReducedWord:
            return False
        return self.main == other.main
    def __hash__(self):
        return repr(self.main).__hash__()
    def draw(self, a = None, size = 70, style = 'til'):
        if style==None:
            style = 'til'
        w = self.main
        n = self.N
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        cols = rainbow(n-1)
        border0 = [-sin(pi*k/n) for k in range(n+1)]
        border1 = [-cos(pi*k/n) for k in range(n+1)]
        p = list(range(1,n+1))
        for k in w:
            a0,b0,c0 = border0[k-1], border0[k], border0[k+1]
            a1,b1,c1 = border1[k-1], border1[k], border1[k+1]
            d0,d1 = a0-b0+c0, a1-b1+c1
            border0[k],border1[k] = d0,d1
            if style == 'til':
                col = cols[p[k]-p[k-1]-1]
                p[k],p[k-1] = p[k-1],p[k]
            elif style == 'rainbow':
                col = cols[k-1]
            else:
                raise Exception('style not defined')
            ax.fill([a0,b0,c0,d0],[a1,b1,c1,d1],facecolor=col,edgecolor='black',linewidth=size/140)
        if a==None:
            plt.show()

def StandardTableau_to_ReducedWord(T):
    T = [list(t) for t in T]
    n = len(T) + 1
    w = []
    for k in range(n*(n-1)//2,0,-1):
        for i in range(n):
            if T[i][-1] == k:
                w.append(i+1)
                x,y = i,n-i-2
                while True:
                    if x>0 and T[x-1][y]!=None:
                        if y>0 and T[x][y-1]!=None:
                            if T[x-1][y]>T[x][y-1]:
                                T[x][y] = T[x-1][y]
                                x,y = x-1,y
                            else:
                                T[x][y] = T[x][y-1]
                                x,y = x,y-1
                        else:
                            T[x][y] = T[x-1][y]
                            x,y = x-1,y
                    elif y>0 and T[x][y-1]!=None:
                        T[x][y] = T[x][y-1]
                        x,y = x,y-1
                    else:
                        T[x][y] = None
                        break
                break
    return W0ReducedWord(n,w)

"""
n = 20
T = StandardTableaux(list(range(n-1,0,-1))).random_element()
w = StandardTableau_to_ReducedWord(T)
w.draw(size = 20)
"""

# S = {(1,2):4, (1,3):2, (2,3):3}
def RewritingGraph(S,w):
    S2 = {}
    for m in S:
        s1,s2 = '',''
        for k in range(S[m]):
            if k%2==0:
                s1 += str(m[0])
                s2 += str(m[1])
            else:
                s1 += str(m[1])
                s2 += str(m[0])
        S2[m] = (s1,s2)
    L = [w]
    E = []
    LM = set()
    while len(L)>0:
        print(len(L))
        L2 = []
        for x in L:
            LM.add(x)
            for i in range(len(x)-1):
                a,b = int(x[i]),int(x[i+1])
                if a<b and x[i:i+S[(a,b)]] == S2[(a,b)][0]:
                    x2 = x[:i] + S2[(a,b)][1] + x[i+S[(a,b)]:]
                    if x2 not in LM:
                        L2.append(x2)
                    E.append((x,x2))
                elif a>b and x[i:i+S[(b,a)]] == S2[(b,a)][1]:
                    x2 = x[:i] + S2[(b,a)][0] + x[i+S[(b,a)]:]
                    if x2 not in LM:
                        L2.append(x2)
                    E.append((x2,x))
                elif a==b:
                    raise Exception("w n'est pas un mot réduit")
        L = L2
    return DiGraph(E)

###### pas abouti : réalisation du graphe des classes de commutations comme somme de Minkowski ####

def VectorProd(x1,x2,y1,y2):
    x3 = sqrt(x1*y1+x2*y2)
    a,b,c = -x3*(x2+y2),x3*(x1+y1),x1*y2-x2*y1
    norm = sqrt(a^2+b^2+c^2)
    return Rational(a),Rational(b),Rational(c)

def Vectors10():
    lv = []
    for k in range(5):
        x1 = cos(pi*(2*k)/5)-cos(pi*(2*k-1)/5)
        x2 = sin(pi*(2*k)/5)-sin(pi*(2*k-1)/5)
        y1 = cos(pi*(2*k+1)/5)-cos(pi*(2*k)/5)
        y2 = sin(pi*(2*k+1)/5)-sin(pi*(2*k)/5)
        lv.append(VectorProd(x1,x2,y1,y2))
        y1 = cos(pi*(2*k+2)/5)-cos(pi*(2*k+1)/5)
        y2 = sin(pi*(2*k+2)/5)-sin(pi*(2*k+1)/5)
        lv.append(VectorProd(x1,x2,y1,y2))
    return lv

def Polytope10():
    lv = Vectors10()
    P = Polyhedron(vertices = [(0,0,0),lv[0]])
    for v in lv[1:]:
        P = P + Polyhedron(vertices = [(0,0,0),v])
    return P

from sage.plot.plot3d.transform import Transformation
from sage.misc.viewer import viewer
viewer.browser('open -a "/Applications/Google Chrome.app"')
viewer.pdf_viewer('evince')

def VectorsView():
    res = point((0, 0, 0), size=0, axes=False, frame=False, projection='orthographic')
    for v in Vectors10():
        res += line( [(0,0,0), v], color="blue", thickness = 10)
    return res


####### classes de symétrie ######

"""
Mots réduits de w0: 			tableaux standards:			pavages:
complémentation + symétrie <->		inv. de schützenberger			symétrie centrale
complémentation <->			conjugué				symétrie d'axe horizontal

pas de classes invariantes par symétrie centrale à partir de A_2
"""

def CommutationMin(w):
    l = w.copy()
    for k in range(1,len(l)):
        t = k
        for i in range(k):
            if l[k-i-1]-1 <= l[k]:
                t = i
                break
        l = l[:k-t] + [l[k]] + l[k-t:k] + l[k+1:]
    return l

def VSCommutationClasses(n):
    W = CoxeterGroup(['A',n-1])
    return [w for w in W.w0.reduced_words() if w==CommutationMin(w)==CommutationMin(list(reversed(w)))]

"""
1, 2, 2, 10, 14 (cf. A180605)
"""


def Alternating_tableau(n):
    T, s = [(n-k)*[0] for k in range(n)], 1
    for k in range(n):
        if k%2==0:
            for i in range(n-k):
                T[k//2][k//2+i] = s
                s += 1
        else:
            for i in range(n-k):
                T[k//2+i+1][k//2] = s
                s += 1
    return StandardTableau(T)