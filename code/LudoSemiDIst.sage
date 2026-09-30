#!/usr/bin/env python
# coding: utf-8

# In[1]:


## Code by Ludovic Schwob (2025)

#def Graph_linear_extension(G):

def Graph_IncidenceMatrix(G):
    n = len(G)
    greater_than = {x:{x} for x in G}
    for x in G:
        gt = greater_than[x]
        for y in G.neighbor_out_iterator(x):
            gt = gt.union(greater_than[y])
        greater_than[x] = gt
    return greater_than

def Graph_IncidenceMatrix(G):
    n = len(G)
    greater_than = {x:{x} for x in G}
    for x in G:
        gt = greater_than[x]
        for y in G.neighbor_out_iterator(x):
            gt = gt.union(greater_than[y])
        greater_than[x] = gt
    return greater_than

def Graph_is_linear_interval(G,x,y):
    it = G.all_paths_iterator([x], [y], simple=True, trivial=True)
    try:
        next(it)
    except StopIteration:  # not comparable
        return False
    try:
        next(it)
    except StopIteration:  # one path
        return True
    return False

def Graph_LinearIntervals(G):
    stock = [(x, x, x) for x in G]
    poly = [len(stock)]
    exposant = 0
    while True:
        exposant += 1
        next_stock = []
        short_stock = [(ch[0], ch[2]) for ch in stock]
        for xmin, cov_xmin, xmax in stock:
            for y in G.neighbor_out_iterator(xmax):
                if exposant == 1:
                    next_stock.append((xmin, y, y))
                elif (cov_xmin, y) in short_stock:
                    if Graph_is_linear_interval(G,xmin, y):
                        next_stock.append((xmin, cov_xmin, y))
        if next_stock:
            poly.append(len(next_stock))
            stock = next_stock
        else:
            break
    return poly


########## for semidistributive lattices ##########

from sage.graphs.distances_all_pairs import distances_all_pairs
def Graph_JoinKMap(G):
    D = distances_all_pairs(G)
    K = {}
    for i in G:
        I = G.neighbors_in(i)
        if len(I)==1:
            for j in I:
                s = False
                for k in G.neighbors_out(j):
                    if D[i][k]==Infinity:
                        s = True
                        I.append(k)
                if not s:
                    break 
            K[i] = j
    return K

def Graph_SDLGraph(G):
    D = distances_all_pairs(G)
    K = {}
    for i in G:
        I = G.neighbors_in(i)
        if len(I)==1:
            for j in I:
                s = False
                for k in G.neighbors_out(j):
                    if D[i][k]==Infinity:
                        s = True
                        I.append(k)
                if not s:
                    break 
            K[i] = j
    return DiGraph({i:[j for j in K if D[i][K[j]]==Infinity] for i in K})

def WeakOrderGraph(W):
    """
    W = ['H',3] par exemple
    """
    W = CoxeterGroup(W)
    L = [W[0]]
    G = {}
    while len(L)>0:
        L2 = set()
        for w in L:
            G[w] = w.upper_covers()
            L2 = L2.union(G[w])
        L = L2
    return DiGraph(G)

def Graph_join_irreducibles_poset(G):
    D = G.distance_all_pairs()
    J = [x for x in G if len(G.neighbors_in(x))==1]
    return Poset((J,lambda p,q : q in D[p]))


######## question Louis ########


#W = CoxeterGroup(W,implementation = 'coxeter3')

degree = lambda w:len(set(w))

"""
Degré sur W tout entier :

type A (A263484, cf. A003319) :
1
1	1
1	2	3
1	3	7	13
1	4	12	32	71

type B (A109281, cf. A109253) :
1	2	5
1	3	9	35
1	4	14	56	309
1	5	20	84	443	3287

type D (A112226, cf. A112225):
1	4	12	40	135
1	5	18	66	293	1537

E_6 :
1	6	25	100	495	3686	47527

F_4 :
1	4	14	76	1057

"""

def ForcingOrderLength(W, size = 20):
    G = WeakOrderGraph(W)
    F = ForcingOrder(Graph_SDLGraph(G))
    F = ForcingOrder(Graph_SDLGraph(G)).relabel(lambda x: tuple(x[0].reduced_word()) )
    D = {}
    for x in F:
        l = degree(x)
        if l in D:
            D[l].append(x)
        else:
            D[l] = [x]
    ll,lc = list(D),rainbow(len(D))
    D = {c: D[l] for l,c in zip(ll,lc)}
    return F.plot(element_colors = D,figsize = size, edge_color = 'lightblue')


# In[32]:


## Code by Ludovic Schwob (2025)

#L = meet-semidistributive lattice
def JoinKMap(L):
    K = {}
    for i in L.join_irreducibles():
        i2 = L.lower_covers(i)[0]
        k = L.join([j for j in L if L.meet(j,i)==i2])
        if L.meet(k,i)!=i2:
            raise Exception("L n'est pas meet-semidistributif")
        K[i] = k
    return K

def SDLGraph(L):
    K = JoinKMap(L)
    return DiGraph({i:[j for j in K if not L.is_lequal(i,K[j])] for i in K})

#G = directed graph, X = subset of G
def right_orthogonal(G,X):
    S = set(X)
    for x in X:
        for y in G.neighbors_out(x):
            S.add(y)
    return [x for x in G if x not in S]

def left_orthogonal(G,X):
    S = set(X)
    for x in X:
        for y in G.neighbors_in(x):
            S.add(y)
    return [x for x in G if x not in S]

#def MOPLattice(G):
#    return LatticePoset(([s for s in Subsets(G) if Set(left_orthogonal(G,right_orthogonal(G,s)))==s],lambda p,q:p.issubset(q)))

#plus performant
def MOPLattice(G,lattice = True):
    LS = [set()]
    L = LS
    while LS!=[]:
        L2 = []
        for s in LS:
            for x in G:
                if x not in s:
                    s2 = s.copy()
                    s2.add(x)
                    t = set(left_orthogonal(G,right_orthogonal(G,s2)))
                    if t not in L and t not in L2:
                        L2.append(t)
        LS = L2
        L.extend(LS)
    if lattice:
        return LatticePoset(([Set(x) for x in L],lambda p,q:p.issubset(q)))
    return L

def OPLattice(G,lattice = True):
    S = set()
    for x in MOPLattice(G,False):
        for s1 in Subsets(Set(x)):
            for s2 in Subsets(Set(right_orthogonal(G,x))):
                S.add((s1,s2))
    if lattice:
        return LatticePoset((S,lambda p,q:p[0].issubset(q[0]) and q[1].issubset(p[1])))
    return S

def SurEdges(G,loops = False):
    E = []
    for x,y,_ in G.edges():
        if loops or x!=y:
            s = True
            for z in G.neighbors_out(y):
                if x!=z and not G.has_edge(x,z):
                    s = False
                    break
            if s:
                E.append((x,y))
    return E

def InEdges(G,loops = False):
    E = []
    for x,y,_ in G.edges():
        if loops or x!=y:
            s = True
            for z in G.neighbors_in(x):
                if z!=y and not G.has_edge(z,y):
                    s = False
                    break
            if s:
                E.append((x,y))
    return E

#expérimental, G = two-acyclic factorization system
#def is_polytopal(G):
#    for x,y,_ in G.edges():
#        if x!=y:
#            s = True
#            for z in G.neighbors_out(y):
#                if not G.has_edge(x,z):
#                    s = False
#                    break
#            if s:
#                for z in G.neighbors_in(x):
#                    if not G.has_edge(z,y):
#                        s = False
#                        break
#                if s:
#                    return False
#    return True

#check if a DiGraph G is a two-acyclic factorization system
def is_TAFS(G):
    GS = DiGraph(SurEdges(G))
    GI = DiGraph(InEdges(G))
    Mult = set((x,y) for x,y,_ in GI.edges())
    for x,y,_ in GS.edges():
        Mult.add((x,y))
        if GS.has_edge(y,x):
#            print(x,y,"cycle d'arêtes surjectives")
            return False
        if y in GI:
            if GI.has_edge(y,x):
#                print(x,y,"cycle d'arêtes surjective et injective")
                return False
            for z in GI.neighbors_out(y):
                Mult.add((x,z))
    for x,y,_ in GI.edges():
        Mult.add((x,y))
        if GI.has_edge(y,x):
#            print(x,y,"cycle d'arêtes injectives")
            return False
    for x,y,_ in G.edges():
        if x!=y and (x,y) not in Mult:
#            print(x,y,'arête composée non présente')
            return False
    return True

def DrawTAFS(G,size = 15):
    S = SurEdges(G)
    I = InEdges(G)
    colI, colS, colIS, col = 'red','blue','black','green'
    D = {colI:[], colS:[], colIS:[], col:[]}
    for e in G.edges():
        if e[:2] in S:
            if e[:2] in I:
                D[colIS].append(e)
            else:
                D[colS].append(e)
        else:
            if e[:2] in I:
                D[colI].append(e)
            else:
                D[col].append(e)
    return G.plot(edge_colors = D,figsize = size)
def DoubleEdgesGraph(G):
    E = []
    for x,y,_ in G.edges():
        if x!=y and G.has_edge(y,x):
            E.append((x,y))
    return Graph(E)

def ForcingOrder(G):
    G = DiGraph(G,loops = True)
    for x in G:
        G.add_edge(x,x)
    S,I = SurEdges(G,True),InEdges(G,True)
    GS,GI = DiGraph(S,loops = True),DiGraph(I,loops = True)
    PS,PI = Poset((G,lambda x,y:(x,y) in S)),Poset((G,lambda x,y:(x,y) in I))
    E = []
    for y in GI:
        for x in PS.subposet(GI.neighbors_in(y)).maximal_elements():
            if x!=y:
                E.append((x,y))
    for y in GS:
        for x in PI.subposet(GS.neighbors_out(y)).minimal_elements():
            if x!=y:
                E.append((x,y))
    D = {x:[] for x in G}
    for x,y in E:
        D[x].append(y)
    D = DiGraph(D)
    F = D.transitive_closure()
    C = {c[0]:c for c in D.strongly_connected_components()}
    for x in C:
        F.merge_vertices(C[x])
    return Poset(F).relabel(lambda x:tuple(C[x]))

def Cov(G,X,PI = None,PS = None):
    if PI==None:
        PI = DiGraph()
        PI.add_vertices(G)
        PI.add_edges(InEdges(G))
        PI = Poset(PI)
    if PS==None:
        PS = DiGraph()
        PS.add_vertices(G)
        PS.add_edges(SurEdges(G))
        PS = Poset(PS)
    return PS.subposet(PI.subposet(X).minimal_elements()).minimal_elements()

def JoinRepresentations(L):
    G = SDLGraph(L)
    PI = DiGraph()
    PI.add_vertices(G)
    PI.add_edges(InEdges(G))
    PI = Poset(PI)
    PS = DiGraph()
    PS.add_vertices(G)
    PS.add_edges(SurEdges(G))
    PS = Poset(PS)
    return {x:Cov(G,[j for j in G if L.is_lequal(j,x)],PI,PS) for x in L}

def SD_Rowmotion(L,lengths = False):
    J = JoinRepresentations(L)
    K = JoinKMap(L)
    A,C = list(L),[]
    while A!=[]:
        a = A[0]
        del A[0]
        c,a = [a],L.meet([K[j] for j in J[a]])
        while a!=c[0]:
            del A[A.index(a)]
            c.append(a)
            a = L.meet([K[j] for j in J[a]])
        C.append(c)
    if lengths:
        L = {}
        for c in C:
            if len(c) in L:
                L[len(c)] += 1
            else:
                L[len(c)] = 1
        return sum(L[k] for k in L),{k:L[k] for k in sorted(L)}
    return C

def WeakOrderJoinRepr(W,w):
    S = W.gens()
    R = []
    for s in w.descents():
        r = w*S[s-1]*w.inverse()
        w2 = w
        while True:
            t = False
            for s2 in w2.descents():
                if w2*S[s2-1]*w2.inverse()!=r:
                    w2 = w2*S[s2-1]
                    t = True
                    break
            if not t:
                R.append(w2)
                break
    return R


def FixedSimpleRoots(W,w):
    S = W.gens()
    return [s for s in S if w*s*w.inverse() in S and (w*s).length()>w.length()]

#pas au point
def WeakIntervalBijection(W,w1,w2):
    wi = w1.inverse()
    wi2 = w2.inverse()
    return (wi,L.join([wi*w2*s*wi2 for s in FixedSimpleRoots(W,w2)]))

def Series1(W):
    R = QQ['x']
    x = R.gen()
    DS = {w:0 for w in W}
    for w in W:
        k = len(FixedSimpleRoots(W,w))
        for w2 in W:
            if w2.weak_le(w):
                DS[w2] += x^k
    return DS

def Series2(W):
    R = QQ['x']
    x = R.gen()
    DS = {w:0 for w in W}
    L = LatticePoset((W,lambda p,q:p.weak_le(q)))
    for w2,w in L.intervals_poset():
        DS[w2.inverse()] += x^Exponent(LatticePoset(L.subposet(L.interval(w2,w))))
    return DS

def LatticeFaces(L,lattice = True):
    LF = []
    for x in L:
        for s in Subsets(L.upper_covers(x)):
            LF.append((x,L.join([x]+list(s))))
    if lattice:
        return Poset((LF,lambda p,q:L.is_lequal(p[0],q[0]) and L.is_lequal(p[1],q[1])))
    return LF

def FaceJoin(L,f1,f2):
    x,y = L.join(f1[0],f2[0]),L.join(f1[1],f2[1])
    l = [y]
    for z in L.lower_covers(y):
        if L.is_lequal(x,z):
            l.append(z)
    return L.meet(l),y

def LatticeFaces2(L):
    G = SDLGraph(L).canonical_label()
    G2 = Derived_Graph(G,[0,2,4,6,8])
    L2 = MOPLattice(G2)
    def face_relabel(X):
        Y = right_orthogonal(G2,X)
        return Set([x for x,k in X if k==0]),Set(left_orthogonal(G,[x for x,k in Y if k==1]))
    return MOPLattice(G),L2.relabel(face_relabel)


#############


def SD_ASMs(n):
    P = Poset((Triples(n),Ordre(Poset({0:[2,3],1:[2,3]}))))
    D = P.hasse_diagram().to_dictionary()
    for x,y in P.intervals_poset():
        if x!=y and not P.covers(x,y):
            if y[0]-x[0]==y[1]-x[1] or y[2]-x[2]==0:
                D[x].append(y)
    G = DiGraph(D)
    for x,y,_ in G.edges():
        if y[0]-x[0]==y[1]-x[1]:
            for z in G.neighbors_in(x):
                if y not in D[z]:
                    D[z].append(y)
        if y[2]-x[2]==0:
            for z in G.neighbors_out(y):
                if z not in D[x]:
                    D[x].append(z)
    return DiGraph(D)     

#G1 et G2 doivent être des relations réflexives
def GraphProd(G1,G2):
    GP = {}
    for x in G1:
        for y in G2:
            lv = []
            for x0 in G1.neighbors_out(x):
                for y0 in G2.neighbors_out(y):
                    lv.append((x0,y0))
            GP[(x,y)] = lv
    return DiGraph(GP)
def GraphSemiDirectProd(G1,G2):
    GP = {}
    for x in G1:
        for y in G2:
            lv = []
            for x0 in G1.neighbors_out(x):
                lv.append((x0,y))
            for y0 in G2.neighbors_out(y):
                if y!=y0:
                    for x0 in G1:
                        lv.append((x0,y0))
            GP[(x,y)] = lv
    return DiGraph(GP)


###############

#l = [0,2,4,8] --> intervalles
#l = [0,2,4,6,8] --> faces
#l = [2,4,8] --> skew intervals

def Derived_Graph(G,l):
    D = {}
    for x in G:
        D[(x,0)] = []
        D[(x,1)] = []
        if 0 in l:
            D[(x,0)].append((x,1))
        if 1 in l:
            D[(x,1)].append((x,0))
        for y in G.neighbors_out(x):
            if x!=y:
                if 2 in l:
                    D[(x,0)].append((y,0))
                if 4 in l:
                    D[(x,0)].append((y,1))
                if 6 in l:
                    D[(x,1)].append((y,0))
                if 8 in l:
                    D[(x,1)].append((y,1))
        for y in G.neighbors_in(x):
            if x!=y:
                if 3 in l:
                    D[(x,0)].append((y,0))
                if 5 in l:
                    D[(x,0)].append((y,1))
                if 7 in l:
                    D[(x,1)].append((y,0))
                if 9 in l:
                    D[(x,1)].append((y,1))
    return DiGraph(D)


def TamariIrreducibles(n):
    I = []
    for k in range(n-1):
        for i in range(1,n-k):
            l = n*[0]
            l[k] = i
            I.append(tuple(l))
    return I

def TamariJoinRepr(v):
    n = len(v)
    l,J = n*[None],[]
    for i,k in enumerate(reversed(v)):
        if k:
            if l[i-k]:
                J.append(tuple((n-i-1)*[0]+[k-l[i-k]]+i*[0]))
            else:
                J.append(tuple((n-i-1)*[0]+[k]+i*[0]))
            l[i-k] = k
    return J

def TamariMeetRepr(v):
    n = len(v)
    l,M = n*[False],[]
    for i,k in enumerate(v):
        if k!=n-i-1 and not l[i+k]:
            M.append(tuple(list(range(n-1,n-i-1,-1))+list(range(k,-1,-1))+list(range(n-i-k-2,-1,-1))))
            l[i+k] = True
    return M

def TamariKMap(v):
    n = len(v)
    K = list(range(n-1,-1,-1))
    for J in TamariJoinRepr(v):
        M = list(range(n-1,-1,-1))
        for i in range(n):
            if J[i]:
                for k in range(i,i+J[i]):
                    M[k] -= n-i-J[i]
                break
        K = [min(i,j) for i,j in zip(K,M)]
    return tuple(K)

def TamariGraph(n):
    K = {i:TamariKMap(i) for i in TamariIrreducibles(n)}
    return DiGraph({i:[j for j in K if any(a>b for a,b in zip(i,K[j]))] for i in K})

def TamariCliques(n):
    G = Graph(TamariGraph(n)).canonical_label()
    G.remove_loops()
    C = [Set(c) for c in sage.graphs.cliquer.all_max_clique(G)]
    return Graph([(C[i],C[j]) for i in range(len(C)) for j in range(i+1,len(C)) if len(C[i].intersection(C[j]))==n-2])

#W = WeylGroup(["A", 4])
#C = list(W.standard_coxeter_elements())
#L = W.cambrian_lattice(C[2]).canonical_label()
#G = SDLGraph(L)
#G.remove_loops()
#Cl = [Set(c) for c in sage.graphs.cliquer.all_max_clique(G)]
#G2 = Graph([(Cl[i],Cl[j]) for i in range(len(Cl)) for j in range(i+1,len(Cl)) if len(Cl[i].intersection(Cl[j]))==3])

#for example if L is a cambrian lattice on A_n, take G = SDLGraph(L)
def SimplicesGraph(G,n):
    G.remove_loops()
    Cl = [Set(c) for c in sage.graphs.cliquer.all_max_clique(G)]
    return Graph([(Cl[i],Cl[j]) for i in range(len(Cl)) for j in range(i+1,len(Cl)) if len(Cl[i].intersection(Cl[j]))==n-1])

#G = SDLGraph(L)
def CliqueSizes(G):
    G.remove_loops()
    Stats(sage.graphs.cliquer.all_cliques(G,1),[lambda c:len(c)])


######

def GraphSkewIntervals(L,lattice = True,l = [2,4,8]):
    G = SDLGraph(L).canonical_label()
    G2 = Derived_Graph(G,l)
    L2 = MOPLattice(G2)
    LL = []
    for x in L2:
        x2 = right_orthogonal(G2,x)
        LL.append((Set(left_orthogonal(G,[a[0] for a in x2 if a[1]==0])),Set([a[0] for a in x if a[1]==1])))
    if lattice:
        return LatticePoset((LL,lambda p,q:p[0].issubset(q[0]) and p[1].issubset(q[1])))
    return LL

def SkewIntervalsInf(L):
    D = {}
    for s in GraphSkewIntervals(L,False):
        if s[0] in D:
            D[s[0]] += 1
        else:
            D[s[0]] = 1
    return D

def SkewIntervalsPivot(L):
    lp = len(L)*[0]
    D = {x:i for i,x in enumerate(L)}
    for s in GraphSkewIntervals(L,False):
        lp[D[s[0].intersection(s[1])]] += 1
    return lp

def SkewIntervalsTest(L,lattice = True):
    G,SI,SF = L.hasse_diagram(),set(),set()
    for x in L:
        V = [L.join(S) for S in Subsets(G.neighbors_out(x))]
        for y in L.order_filter([x]):
            for z in V:
                SI.add((z,y))
        V = [L.meet(S) for S in Subsets(G.neighbors_in(x))]
        for y in L.order_ideal([x]):
            for z in V:
                SF.add((y,z))
    LL = [x for x in SI if x in SF]
    if lattice:
        return LatticePoset((LL,lambda p,q:L.is_lequal(p[0],q[0]) and L.is_lequal(p[1],q[1])))
    return LL

def Exponent(L):
    K,C = JoinKMap(L),set(L.coatoms())
    return sum(1 for c in L.atoms() if K[c] in C)

def SkewIntervals(L,lattice = True):
    SI = []
    for x,y in L.intervals_poset():
        I = LatticePoset(L.subposet(L.interval(x,y)))
        K,C = JoinKMap(I),set(I.coatoms())
        cc = [(c,K[c]) for c in I.atoms() if K[c] in C]
        for S in Subsets(cc):
            SI.append((I.join([s[0] for s in S]),I.meet([s[1] for s in S])))
    if lattice:
        return LatticePoset((SI,lambda p,q:L.is_lequal(p[0],q[0]) and L.is_lequal(p[1],q[1])))
    return SI

def SkewFaces(L,lattice = True):
    SF = []
    for x in L:
        for J in Subsets(L.upper_covers(x)):
            y = L.join([x]+list(J))
            I = LatticePoset(L.subposet(L.interval(x,y)))
            K,C = JoinKMap(I),set(I.coatoms())
            cc = [(c,K[c]) for c in I.atoms() if K[c] in C]
            for S in Subsets(cc):
                SF.append((I.join([s[0] for s in S]),I.meet([s[1] for s in S])))
    if lattice:
        return LatticePoset((SF,lambda p,q:L.is_lequal(p[0],q[0]) and L.is_lequal(p[1],q[1])))
    return SF

#pas intéressant
def NilpotentIntervals(L):
    LI = L.intervals_poset()
    D = {I:0 for I in LI}
    for x,y in LI:
        I = LatticePoset(L.subposet(L.interval(x,y)))
        K,C = JoinKMap(I),set(I.coatoms())
        cc = [c for c in I.atoms() if K[c] in C]
        D[(I.join(cc),y)] += (-1)**len(cc)
    return D

def DPolyIntervals(LI):
    """
    par exemple
    L = posets.TamariLattice(n).canonical_label()
    LI = L.intervals_poset()
    ou
    LI = SkewIntervals(L)
    """
    R4.<x,y,z,w> = ZZ['x,y,z,w']
    D = R4(0)
    for i in LI:
        m = 1
        for j in LI.upper_covers(i):
            if i[0]==j[0]:
                m *= y
            elif i[1]==j[1]:
                m *= x
            else:
                print('erreur1')
        for j in LI.lower_covers(i):
            if i[0]==j[0]:
                m *= w
            elif i[1]==j[1]:
                m *= z
            else:
                print('erreur2')
        D += m
    return D


#Stats(L.intervals_poset(),[lambda I:Exponent(L.sublattice(L.interval(I[0],I[1])))])

#(v1,v2) = deux brackets vectors délimitant un intervalle
def TamariExponent(v1,v2):
    n = len(v1)
    l1 = n*[None]
    for i,k in enumerate(v1):
        if k!=n-i-1 and l1[i+k]!='f':
            l1[i] = k+1
            if k>0:
                l1[i+k] = 'f'
    l2 = list(v2)
    for i,k in enumerate(reversed(v2)):
        if k:
            l2[n-i-1] = k-l2[n-i-1+k]
            l2[n-i-1+k] = k
    return sum(1 for i in range(n) if l1[i]==l2[i] and v2[i]!=0)

def possible_values(bv):
    n = len(bv)
    L,i = [0],0
    while i < n:
        L.append(bv[i]+i+1)
        i += bv[i]+1
    return L 
def TamariIntervals(n):
    L=[]
    if n == 1:
        return([([0],[0])])
    for bvs in TamariIntervals(n-1):
        for b in possible_values(bvs[1]):
            for a in possible_values(bvs[0]):
                if b>=a:
                    L.append(([a]+bvs[0],[b]+bvs[1]))
    return L


######### ANY LATTICE #########@

#M = 0,1 matrix encoding a relation between two sets X Y, and v a vector encoding a subset of X
def right_orthogonal2(M,v):
    return vector([int(k==0) for k in M*v],immutable = True)
#v encods a subset of Y
def left_orthogonal2(M,v):
    return vector([int(k==0) for k in v*M],immutable = True)

def Lattice(M,lattice = True):
    n = len(M[0])
    LS = [vector(n*[0],immutable = True)]
    L = LS
    while LS!=[]:
        L2 = []
        for s in LS:
            for x in range(n):
                if s[x]==0:
                    s2 = list(s)
                    s2[x] = 1
                    s2 = vector(s2,immutable = True)
                    t = left_orthogonal2(M,right_orthogonal2(M,s2))
                    if t not in L and t not in L2:
                        L2.append(t)
        LS = L2
        L.extend(LS)
    if lattice:
        return LatticePoset((L,lambda p,q:all(i<=j for i,j in zip(p,q))))
    return L

####### tests ######

def RandomOrientation(G):
    G.remove_loops()
    E,E2 = G.edges(),[]
    for e in E:
        if randint(0,1):
            E2.append(e)
        else:
            E2.append((e[1],e[0]))
    return DiGraph(E2)

#E = liste d'arêtes
def RandomOrientations(E,n = 1000,LL = set()):
    for k in range(n):
        E2 = []
        for e in E:
            if randint(0,1):
                E2.append(e)
            else:
                E2.append((e[1],e[0]))
        G = DiGraph(E2)
        if is_TAFS(G):
            LL.add(G.canonical_label().copy(immutable = True))
    return LL

"""
Plus petit treillis semidistributif mais pas congruence-uniforme:
L = LatticePoset({0:[1,2],1:[3,4],2:[4,5],3:[6,7],4:[8],5:[9],6:[11],7:[8,10],8:[9],9:[12],10:[11,12],11:[13],12:[13]})

"""

