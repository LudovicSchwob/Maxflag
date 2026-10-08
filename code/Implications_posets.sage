def MaxflagCompatibilities(n):
    L = LatticePoset(minMaxFlag(n).poset())
    G = Graph(is_maxflag(L, return_complex = True)[1])
    l = []
    for p, q, _ in G.edges():
        a = p.nonCrossingArcDiag().arc_list[0]
        b = q.nonCrossingArcDiag().arc_list[0]
        l.append(ArcDiag([a, b], n))
    return l

def arc_to_triplet(a, translate = False):
    if translate:
        n = a.n
        return a.i - (n+2)/4, a.i+len(a.L) - n/2, a.j-1 - (3*n-2)/4
    return a.i, a.i+len(a.L), a.j-1

from sage.plot.plot3d.transform import Transformation
from sage.misc.viewer import viewer
viewer.browser('open -a "/Applications/Google Chrome.app"')
viewer.pdf_viewer('evince')


def Maxflag3D(n):
    L = LatticePoset(minMaxFlag(n).poset())
    G = Graph(is_maxflag(L, return_complex = True)[1])
    res = point((0, 0, 0), size=0, axes=False, frame=False, projection='orthographic')
    vertices, edges = [], []
    J = L.join_irreducibles()
    for i, p in enumerate(J):
        t = arc_to_triplet(p.nonCrossingArcDiag().arc_list[0], True)
        vertices.append(t)
        res += point(t, color = "black", size=10)
        for q in J[i+1:]:
            if G.has_edge(p, q):
                u = arc_to_triplet(q.nonCrossingArcDiag().arc_list[0], True)
                res += line([t, u], color= "blue", thickness = 2)
    return res.transform(T=Transformation(m=matrix(4,4,[[0,-1,0,0],[-1,1,-1,0],[1,0,-1,0],[0,0,0,1]])))

################################################################################################

def PosetUnionDif(P1, P2, Q):
    S1, S2 = set([tuple(e) for e in P1.relations()]), set([tuple(e) for e in P2.relations()])
    T = set([tuple(e) for e in Q.relations()])
    U = []
    for e in S1:
        if e in S2:
            if e not in T:
                print('Arête dans P1 et P2 mais pas dans Q !')
            U.append(e)
        if e not in T:
            U.append(e)
    for e in S2:
        if e not in S1 and e not in T:
            U.append(e)
    G = DiGraph(U, loops = True)
    G.remove_loops()
    return Poset(G)

def RectanglePart(p, i0 = 1):
    if len(p)==0:
        return []
    return RectanglePart([k for k in p if k<p[0]],i0) + [[(i,j) for i in range(i0, p[0]+1) for j in range(p[0]+1, len(p)+i0+1)]] + RectanglePart([k for k in p if k>p[0]],p[0]+1)

def RectanglesPoset(p):
    P = Poset()
    for l in RectanglePart(p):
        P = P.disjoint_union(Poset((l,lambda i,j:i[0]>=j[0] and i[1]<=j[1]))).relabel(lambda x:x[1])
    return P

def JuraCambPosets(n):
    return [RectanglesPoset([k] + list(range(k-1, 0, -1)) + list(range(k+1, n))) for k in range(1, n)]

def inversionSet(p):
    return [(j, i) for k, i in enumerate(p) for j in p[k+1:] if i > j]

def IdealInversions(P, n):
    I = {}
    for p in Permutations(n):
        S = Set(inversion_set(p))
        if P.is_order_ideal(S):
            I[p] = S
    return LatticePoset((I, lambda p,q: I[p].issubset(I[q])))


def RootPoset(n):
    return Poset(([(i, j) for i in range(1, n) for j in range(i+1, n+1)], lambda p, q: p[0] >= q[0] and p[1] <= q[1]))

def ImplicationPoset(Q):
    """
    Q = object of type WOQuotient
    This construction only makes sense for permutree quotients
    For other quotients (with arcs L,i,j,R if j-i > 2),
    we would need implications like "x AND y imply z" which is shit
    """
    n = Q.n
    D = {Permutation(p.perm): set([(i,j) for j,i in p.invSet()]) for p in Q.minElements(perm_object = True)}
    Phi = RootPoset(n)
    G = DiGraph()
    for x in Phi:
        G.add_vertex(x)
    for x, y in Phi.relations():
        if x != y:
            if all(y not in s or x in s for s in D.values()):
                G.add_edge(x, y)
    P = Poset(G)
    for p in Permutations(n):
        s = set(inversionSet(p))
        if P.is_order_ideal(s) and not p in D:
            print('Pas dans le quotient:', p)
    return P

def PermutreePoset(d):
    """
    d = decorations of a permutree coded as integers 0...3
    """
    n = len(d) + 2
    G = DiGraph()
    for i in range(1, n):
        for j in range(i+1, n+1):
            G.add_vertex((i, j))
    for k, c in enumerate(d):
        if c%2:
            for i in range(1, k+2):
                for j in range(k+3, n+1):
                    G.add_edge(((k+2, j), (i, j)))
        if c//2:
            for i in range(1, k+2):
                for j in range(k+3, n+1):
                    G.add_edge(((i, k+2), (i, j)))
    return Poset(G)

def left_packing(X, d, up = True):
    """
    X = subset of root poset (as a list of pairs (i, j) with i < j)
    d = decorations of a permutree coded as integers 0...3
    """
    n = len(d) + 2
    S = set(X)
    for k, c in enumerate(d):
        if c//2:
            for i in range(1, k+2):
                x = (i, k+2)
                if (up and x in S) or (not up and x not in S):
                    for j in range(k+3, n+1):
                        y = (i, j)
                        if up and y not in S:
                            S.remove(x)
                            S.add(y)
                            break
                        elif not up and y in S:
                            S.remove(y)
                            S.add(x)
                            break
    return S

def right_packing(X, d, up = True):
    """
    X = subset of root poset (as a list of pairs (i, j) with i < j)
    d = decorations of a permutree coded as integers 0...3
    """
    n = len(d) + 2
    S = set(X)
    for k in range(n-3,-1,-1):
        c = d[k]
        if c%2:
           for j in range(k+3, n+1):
                x = (k+2, j)
                if (up and x in S) or (not up and x not in S):
                    for i in range(k+1, 0 ,-1):
                        y = (i, j)
                        if up and y not in S:
                            S.remove(x)
                            S.add(y)
                            break
                        elif not up and y in S:
                            S.remove(y)
                            S.add(x)
                            break
    return S

"""
Permutree lattices satisfying this condition are those such that in between each pair of decorations 3 (X), there is at least a 0 (|).
Their number is given by
1, 4, 15, 55, 200, 725... (OEIS A039717)
Satisfies the following recursion
(by distinguishing according to "a_n : a 0 appears after the last 3", and "b_n : the contrary") 
a_0 = 1, b_1 = 0, a_{n+1} = 3*a_n + b_n, b_{n+1} = a_n + 2*b_n
"""
def testLeftRight(d, up = True, counter_example = False):
    n = len(d) + 2
    for p in Permutations(n):
        X = inversionSet(p)
        if left_packing(right_packing(X, d, up), d, up) != right_packing(left_packing(X, d, up), d, up):
            if counter_example:
                print(p, X)
            return False
    return True

def packingInvariants(n):
    A, B = [[]], []
    for _ in range(n-2):
        A, B = [m + [k] for m in A for k in range(3)] + [m + [0] for m in B], [m + [3] for m in A] + [m + [k] for m in B for k in range(1,3)]
    return A + B
    
# injectif si et seulement si les décorations sont dans [0, 1, 2]
# quid des autres treillis ???? tailles ????
def testInjectivePerm(n, up = True):
    l = []
    for d in packingInvariants(n):
        s, t = set(), True
        for p in Permutations(n):
            S = Set(right_packing(left_packing(inversionSet(p), d, up), d, up))
            if S in s:
                t = False
                break
            s.add(S)
        if t:
            l.append(d)
    return l

"""
Le treillis obtenu semble être meet-distributif si et seulement si d \in [1,2]^n
 Dans ce cas là, est-ce le dual d'un sorting order ????? (apparement oui)
 - join-semidistributif :  d \in [1,2]^n  OU  d est de la forme 1^i 0^j 2^k
 - semidistributif : si d est de la forme 1^i 0^j 2^k
 - meet-semidistributif si d n'a pas de 1 après un 2
 - join-distributif / distributif si de la forme 1^i 2^j
"""
def PermutreeMiddleOrder(d, up = True):
    """
    d = tuple à valeurs dans [0,1,2]
    """
    n = len(d) + 2
    D = {p: Set(right_packing(left_packing(inversionSet(p), d, up), d, up)) for p in Permutations(n)}
    return LatticePoset((D, lambda p,q: D[p].issubset(D[q])))

# c'est bien des middle orders, mais pas toujours semidistributifs
def testMiddleOrders(n, up = True):
    semi = []
    for d in Tuples([0,1,2],n-2):
        L = PermutreeMiddleOrder(d, up)
        if L.is_join_distributive():
            semi.append(d)
#        for p, q in L.cover_relations():
#            if not p.bruhat_lequal(q):
#                print('pas inclus dans Bruhat', d, p, q)
#                break
    return semi

"""
Les middle orders max-flag semblent être ceux de la forme
1^i2^j ou 1^i02^j (i,j >=0)
Les seuls semidistributifs maxflags sont alors les jurassiens cambriens

Les middle orders dont le dual est max-flag semblent être ceux de la forme
2^i1^j, 1^i2^j, ou 1^i2^j01^i2^j (i,j >=0)

	8, 16, 28, 45, 68, 98...(A255993)
dont	4, 6,  8,  10, 12, 14...  semidistributifs
"""
def MaxflagMiddleOrders(n, up = True, dual = False):
    l = []
    for d in Tuples([0,1,2],n-2):
        L = PermutreeMiddleOrder(d, up)
        if dual:
            L = L.dual()
        if is_maxflag(L):
            l.append(d)
    return l

def testSortingOrders(n):
    W = WeylGroup(['A', n-1], prefix = 's')
    D = {tuple(w): WSortingOrder(W, w).relabel(lambda x: Permutation(list(x.to_permutation()))) for w in CommutationClasses(W)}
    S = {}
    for d in Tuples([1, 2], n-2):
        L = PermutreeMiddleOrder(d)
        for w in D:
            if L == D[w]:
                print(d, w)
                S[d] = W0ReducedWord(n, w)
    return S

"""
Sorting words corresponding to cambrian middle orders :
(1, 1, 1) (1, 2, 1, 3, 2, 1, 4, 3, 2, 1)
(2, 1, 1) (2, 1, 3, 2, 1, 4, 3, 2, 1, 4)
(1, 2, 1) (1, 3, 2, 1, 4, 3, 2, 1, 4, 3)
(2, 2, 1) (3, 2, 1, 4, 3, 2, 1, 4, 3, 4)
(1, 1, 2) (1, 2, 1, 4, 3, 2, 1, 4, 3, 2)
(2, 1, 2) (2, 1, 4, 3, 2, 1, 4, 3, 2, 4)
(1, 2, 2) (1, 4, 3, 2, 1, 4, 3, 2, 4, 3)
(2, 2, 2) (4, 3, 2, 1, 4, 3, 2, 4, 3, 4)
"""

########### generalisation autres types #############

def cambrian_min(W, w, c):
    if len(c)==0 or w==W(1):
        return w
    S = W.gens()
    s = S[c[0]-1]
    if (s*w).length() < w.length():
        return s * cambrian_min(W, s*w, c[1:]+[c[0]])
    c = c[1:]
    p = W(1)
    while True:
        t = False
        for k in w.descents('left'):
            if k in c:
                p *= S[k-1]
                w = S[k-1]*w
                t = True
                break
        if not t:
            return cambrian_min(W, p, c)

def cambrian_max(W, w, c):
    return cambrian_min(W, w*W.w0, list(reversed(c)))*W.w0

def Cambrian_Quotient(W, c, maxima = True):
    """
    c = mot de coxeter, par exemple:
    W = WeylGroup(['D',4],prefix = 's')
    c = [3,1,4,2]
    Cambrian_Quotient(W, c)
    """
    if maxima:
        return {w: cambrian_max(W, w, c) for w in W}
    return {w: cambrian_min(W, w, c) for w in W}



def Cambrian_Lattice(W, c, maxima = True):
    if maxima:
        return LatticePoset(([w for w in W if w == cambrian_max(W, w, c)], lambda p,q: p.weak_le(q)))
    return LatticePoset(([w for w in W if w == cambrian_min(W, w, c)], lambda p,q: p.weak_le(q)))

"""
IL FAUDRAIT PEUT-ÊTRE FAIRE UN TASSAGE VERS LE HAUT POUR OBTENIR UN MIDDLE ORDER QUI EST UN SORTING ORDER

"""


# FAIRE LA MÊME CHOSE POUR LES PERMUTARBRES OÙ UN SEUL ARC EST CONTRACTÉ
# TROUVER SI ÇA CORRESPOND À DES MIDDLE ORDERS
def ImplicationPoset(W, c, maxima = True):
    if type(W) == list:
        W = WeylGroup(W, prefix = 's')
    Phi = RootSystem(W).root_poset()
#    Phi = Poset(([x.element for x in Phi], lambda p, q: Phi.is_lequal(p, q)))
    G = DiGraph()
    for x in Phi:
        G.add_vertex(x.element)
    L = Cambrian_Lattice(W, c, maxima)
    D = {w: set(w.inversions(side = 'left', inversion_type = 'roots')) for w in L}
    for x, y in Phi.relations():
        if x != y:
            if maxima:
                if all(x not in s or y in s for s in D.values()):
                    G.add_edge(x, y)
            else:
                if all(y not in s or x in s for s in D.values()):
                    G.add_edge(x, y)
    P = Poset(G)
    for w in W:
        s = w.inversions(side = 'left', inversion_type = 'roots')
        if maxima:
            if P.is_order_filter(s) and not w in D:
                print('Pas dans le quotient:', w)
        else:
            if P.is_order_ideal(s) and not w in D:
                print('Pas dans le quotient:', w)
    return P