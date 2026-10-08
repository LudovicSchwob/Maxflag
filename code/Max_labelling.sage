# teste si pour un treillis L, pour tout x dans L |cjr(x)] <= |mjr(x)|
def TestCjrMjr(L):
    P = L.join_irreducibles_poset()
    for x in L:
        m = P.subposet([j for j in P if L.is_lequal(j,x)]).maximal_elements()
        if len(m) < len(L.canonical_joinands(x)):
            print(x)
            return False
    return True

def Maxlabelling(L):
    D = {}
    P = L.join_irreducibles_poset()
    for y in L:
        m = P.subposet([j for j in P if L.is_lequal(j,y)]).maximal_elements()
        for x in L.lower_covers(y):
            i = Edge_to_JoinIrr(L, x, y)
            for j in m:
                if P.is_lequal(i, j):
                    D[(x, y)] = j
    return D

"""
Les max-congruences de Tamari sont comptées par Fibonacci_{2n} !
"""

def MaxForcingOrder(L):
    D = Maxlabelling(L)
    F = set()
    for y in L:
        for s in Subsets(L.lower_covers(y), 2):
            a, b = tuple(s)
            d1, d2 = D[(a, y)], D[(b, y)]
            x = L.meet([a, b])
            I = L.subposet(L.interval(x,y))
            for c, d in I.cover_relations():
                if d not in [a, b, y]:
                    if D[(c, d)] != d1:
                        F.add((d1, D[(c, d)]))
                    if D[(c, d)] != d2:
                        F.add((d2, D[(c, d)]))
    G = DiGraph()
    for x in L.join_irreducibles():
        G.add_vertex(x)
    for x, y in F:
        G.add_edge(x, y)
    return G

def MaxQuotient(L, contracted_irr):
    D = Maxlabelling(L)
    return LatticePoset(L.subposet([y for y in L if all(D[(x,y)] not in contracted_irr for x in L.lower_covers(y))]))

"""
Quelques treillis semblent avoir tous leurs max-quotients semidistributifs :
 - Les jurassiens dont la rectangulation est une union de chaînes
 - les Alt-Tamari avec delta = 1^i0^j
"""

def TestQuotient(L):
    F = Poset(MaxForcingOrder(L))
    for s in F.antichains():
        try:
            Q = MaxQuotient(L, F.order_filter(s))
#            if not Q.is_join_semidistributive():
            if not is_maxflag(Q):
                print(s)
        except:
            print('bug', s)