#ref: Drew Armstrong, The sorting order on a Coxeter group

# A TESTER : sorting orders sur les groupes de Coxeter affine
# tester qqch entre SortingOrder et InfiniteSortingOrder, en prenant
# w = w_J * ^Jw^\infty où w_J est un sorting word pour W_J, et ^Jw un mot de Coxeter par exemple
# notamment dans le cas où W_J est le groupe de Weyl correspondant

#par exemple W = WeylGroup(['A', 3],prefix = 's')
#par ex. w = [1,3,2,1,2]
#x = élément de W
def SortingSubword(W,w,x):
    s,y = [],x
    for i,k in enumerate(w):
        if y.has_left_descent(k):
            s.append(i)
            y = y.apply_simple_reflection_left(k)
    return s

#w = mot de coxeter (par ex. [3,1,2], contient une fois chaque générateur)
#sorting word = wwwwwwww...
def InfiniteSortingSubword(W,w,x):
    s,y,j = [],x,0
    while y.length()>0:
        for i,k in enumerate(w):
            if y.has_left_descent(k):
                s.append(i+len(w)*j)
                y = y.apply_simple_reflection_left(k)
        j += 1
    return s

def SortingOrder(W,w):
    return LatticePoset(([Set(SortingSubword(W,w,x)) for x in W],lambda p,q:p.issubset(q)))
def InfiniteSortingOrder(W,w):
    return LatticePoset(([Set(InfiniteSortingSubword(W,w,x)) for x in W],lambda p,q:p.issubset(q)))

def CSortableOrder(W,w):
    return LatticePoset(([Set(InfiniteSortingSubword(W,w,x)) for x in W if x.is_coxeter_sortable(w)],lambda p,q:p.issubset(q)))

#w = reduced word for w0 to get an order on W
#reduced words give the same lattice if and only if they are in the same commutation class
def WSortingOrder(W,w):
    S = W.gens()
    def relabel(x):
        r = W(1)
        for k in sorted(x):
            r *= S[w[k]-1]
        return r
    return SortingOrder(W,w).relabel(relabel)

def CommutationClasses(W):
    G = W.coxeter_diagram()
    return [w for w in W.w0.reduced_words() if all(w[i-1]<w[i] or G.has_edge(w[i-1],w[i]) for i in range(1,len(w)))]

def CommutationMin(W, w):
    G = W.coxeter_diagram()
    w = w.copy()
    k, n = 0, len(w)
    while True:
        if k!=0 and w[k-1]>w[k] and not G.has_edge(w[k-1],w[k]):
            w[k-1], w[k] = w[k], w[k-1]
            k -= 1
        elif k!=n-1:
            if w[k]>w[k+1] and not G.has_edge(w[k],w[k+1]):
                w[k], w[k+1] = w[k+1], w[k]
            k += 1
        else:
            return w

def WSortingOrders(W):
    LL = {}
    for w in CommutationClasses(W):
        L = WSortingOrder(W,w)
        if L not in LL:
            LL[L] = w
    return LL

################ sorting orders distributifs

"""
en type B: 2, 6, 19

[2, 1, 2, 3, 2, 1, 3, 2, 3]
[1, 2, 1, 3, 2, 1, 3, 2, 3]

[1, 3, 2, 3, 1, 2, 1, 3, 2]
[1, 3, 2, 3, 2, 1, 2, 3, 2]

[3, 2, 3, 2, 1, 2, 3, 2, 1]
[2, 3, 2, 3, 1, 2, 3, 2, 1]


[1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 4, 3, 4]
[1, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 4, 3, 4]
[2, 3, 2, 1, 2, 3, 4, 3, 2, 1, 4, 3, 2, 4, 3, 4]
[3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 4, 3, 2, 4, 3, 4]
[2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 4, 3, 4]

[1, 2, 1, 4, 3, 4, 3, 2, 3, 4, 3, 1, 2, 3, 4, 3]
[1, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 4, 3]
[2, 1, 2, 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 4, 3]
[2, 1, 2, 4, 3, 4, 3, 2, 1, 3, 4, 3, 2, 3, 4, 3]

[1, 4, 3, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2]
[1, 4, 3, 4, 3, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2]
[1, 3, 4, 3, 4, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2]
[1, 3, 4, 3, 4, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2]

[3, 4, 3, 2, 4, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1]
[2, 4, 3, 4, 3, 2, 3, 4, 3, 1, 2, 3, 4, 3, 2, 1]
[4, 3, 4, 3, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1]
[2, 4, 3, 2, 4, 3, 2, 4, 3, 1, 2, 3, 4, 3, 2, 1]
[2, 3, 2, 4, 3, 2, 4, 3, 4, 1, 2, 3, 4, 3, 2, 1]
[3, 2, 3, 4, 3, 2, 4, 3, 4, 1, 2, 3, 4, 3, 2, 1]

en type D:

[2, 1, 2, 3, 2, 1, 4, 2, 1, 3, 2, 4]
[2, 3, 2, 1, 2, 3, 4, 2, 1, 3, 2, 4]
[3, 2, 1, 3, 2, 3, 4, 2, 1, 3, 2, 4]
[1, 2, 1, 3, 2, 1, 4, 2, 1, 3, 2, 4]
[1, 3, 2, 1, 3, 2, 4, 2, 1, 3, 2, 4]

[2, 4, 2, 1, 2, 4, 3, 2, 1, 4, 2, 3]
[2, 1, 2, 4, 2, 1, 3, 2, 1, 4, 2, 3]
[1, 2, 1, 4, 2, 1, 3, 2, 1, 4, 2, 3]
[4, 2, 1, 4, 2, 4, 3, 2, 1, 4, 2, 3]
[1, 4, 2, 1, 4, 2, 3, 2, 1, 4, 2, 3]

[2, 4, 2, 3, 2, 4, 1, 2, 3, 4, 2, 1]
[3, 4, 2, 3, 4, 2, 1, 2, 3, 4, 2, 1]
[3, 2, 3, 4, 2, 3, 1, 2, 3, 4, 2, 1]
[4, 2, 3, 4, 2, 4, 1, 2, 3, 4, 2, 1]
[2, 3, 2, 4, 2, 3, 1, 2, 3, 4, 2, 1]

Il semblerait que les seuls middle orders en type D soient les middle orders minuscules !
Testé jusqu'à D_5

F_4 :
[1, 3, 2, 1, 3, 2, 1, 3, 2, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4]
[2, 1, 2, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4]
[2, 3, 2, 1, 3, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4]
[1, 3, 2, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4]
[1, 2, 1, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4]
[2, 1, 2, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4]
[2, 3, 2, 1, 3, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4]
[1, 3, 2, 1, 3, 2, 1, 3, 2, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4]
[1, 2, 1, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4]
[1, 3, 2, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4]
[3, 2, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4]
[3, 2, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4]

[2, 4, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1]
[2, 4, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1]
[4, 3, 2, 4, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]
[3, 2, 3, 2, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1]
[2, 3, 2, 3, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1]
[4, 3, 2, 4, 3, 2, 4, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1]
[2, 4, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]
[2, 4, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]
[3, 2, 3, 2, 4, 3, 2, 1, 3, 2, 4, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]
[2, 3, 2, 3, 4, 3, 2, 1, 3, 2, 4, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]
[3, 4, 3, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1]
[3, 4, 3, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]

H_3:

[3, 2, 3, 2, 3, 1, 2, 3, 2, 1, 3, 2, 3, 2, 1]
[2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2, 3, 2, 1]

[2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2, 3]
[2, 1, 2, 3, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 3]
[1, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 3]
[1, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2, 3]
[2, 1, 2, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 3]
[1, 2, 1, 3, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 3]

[1, 3, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2]
[1, 3, 2, 3, 2, 1, 3, 2, 3, 2, 1, 3, 2, 3, 2]

H_4: Gott weiß
"""

def ParabolicQuotient(W,J):
    try:
        W = WeylGroup(W,prefix = 's')
    except:
        W = CoxeterGroup(W)
    L = [W[0]]
    l = L
    while len(l)>0:
        l2 = []
        for x in l:
            for y in x.upper_covers():
                if any(k not in J for k in (x*y.inverse()).reduced_word()) and y not in l2:
                    l2.append(y)
        l = l2
        L.extend(l2)
    return L

def ParabolicQuotientSortingWords(W, J, omega):
    """
    omega must be a reduced word for w0 of W_J yielding a distributive sorting order
    """
    Q = ParabolicQuotient(W, J)
    W = CoxeterGroup(W)
    G = W.coxeter_diagram()
    lw = []
    for w in Q[-1].reduced_words():
        if all(w[i-1]<w[i] or G.has_edge(w[i-1],w[i]) for i in range(1,len(w))):
            L = SortingOrder(W, omega + w)
            if L.is_distributive():
                lw.append(w)
    return lw

@cached_function
def A_SortingWords(n):
    if n==0:
        return [[]]
    l = []
    for k in range(n):
        w = []
        for i in range(n-k):
            w.extend(range(i+k+1,i,-1))
        for w1 in A_SortingWords(k):
            for w2 in A_SortingWords(n-k-1):
                l.append(w1 + [i+k+1 for i in w2] + w)
    return l

@cached_function
def B_SortingWords(n):
    if n==0:
        return [[]]
    l = []
    for k in range(n):
        w = []
        for i in range(k+1):
            w.extend(range(k+1,n))
            w.extend(range(n,i,-1))
        W = [w]
        if k!=0 and k!=n-1:
            w = []
            for i in range(1,n-k):
                w.extend(range(k+i,i-1,-1))
            for i in range(k+1):
                w.extend(range(n,i,-1))
            W.append(w)
        for w1 in A_SortingWords(k):
            for w2 in B_SortingWords(n-k-1):
                for w in W:
                    l.append(w1 + [i+k+1 for i in w2] + w)
    return l

@cached_function
def B_SortingWords2(n):
    if n==2:
        return [[1,2,1,2],[2,1,2,1]], []
    l1,l2 = [], []
    for k in range(n):
        w = []
        for i in range(k+1):
            w.extend(range(k+1,n))
            w.extend(range(n,i,-1))
        W = [w]
        if k!=0 and k!=n-1:
            w = []
            for i in range(1,n-k):
                w.extend(range(k+i,i-1,-1))
            for i in range(k+1):
                w.extend(range(n,i,-1))
            W.append(w)
        for w1 in A_SortingWords(k):
            if k<n-2:
                m1, m2 = B_SortingWords2(n-k-1)
                for w2 in m1:
                    for w in W:
                        l1.append(w1 + [i+k+1 for i in w2] + w)
                for w2 in m2:
                    for w in W:
                        l2.append(w1 + [i+k+1 for i in w2] + w)
            else:
                for w2 in B_SortingWords(n-k-1):
                    for w in W:
                        l2.append(w1 + [i+k+1 for i in w2] + w)
    return l1, l2

@cached_function
def D_SortingWords(n):
    if n==3:
        relabel = {1: 2, 2: 1, 3: 3}
        return [[relabel[k] for k in w] for w in A_SortingWords(3)]
    l = []
    w = []
    for i in range(n-1):
        if i%2==0:
            w.append(n)
        else:
            w.append(n-1)
        w.extend(range(n-2,i,-1))
    for w1 in A_SortingWords(n-1):
        l.append(w1 + w)
        l.append([i if i<n-1 else 2*n-i-1 for i in w1 + w])
    w = list(range(1,n+1)) + list(range(n-2,0,-1))
    for w1 in D_SortingWords(n-1):
        l.append([i+1 for i in w1] + w)
    return l

def F4_SortingWords():
    l1 = [[2, 3, 2, 1, 3, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4],
          [2, 3, 2, 1, 3, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4],
          [3, 2, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4],
          [3, 2, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4],
          [3, 2, 3, 2, 4, 3, 2, 3, 4, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1],
          [2, 3, 2, 3, 4, 3, 2, 3, 4, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1],
          [3, 2, 3, 2, 4, 3, 2, 3, 4, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1],
          [2, 3, 2, 3, 4, 3, 2, 3, 4, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]]
    l2 = [[1, 3, 2, 1, 3, 2, 1, 3, 2, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4],
          [2, 1, 2, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4],
          [1, 3, 2, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4],
          [1, 2, 1, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 3, 4],
          [2, 1, 2, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4],
          [1, 3, 2, 1, 3, 2, 1, 3, 2, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4],
          [1, 2, 1, 3, 2, 1, 3, 2, 3, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4],
          [1, 3, 2, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1, 4, 3, 2, 3, 4],
          [2, 4, 3, 2, 4, 3, 2, 4, 3, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1],
          [2, 4, 3, 2, 3, 4, 3, 2, 3, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1],
          [4, 3, 2, 4, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1],
          [4, 3, 2, 4, 3, 2, 4, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1],
          [2, 4, 3, 2, 3, 4, 3, 2, 3, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1],
          [2, 4, 3, 2, 4, 3, 2, 4, 3, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1],
          [3, 4, 3, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 1, 4, 3, 2, 1, 3, 2, 4, 3, 2, 1],
          [3, 4, 3, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1, 2, 3, 2, 4, 3, 2, 1]]
    return l1, l2

def H3_SortingWords():
    l1 = [[3, 2, 3, 2, 3, 1, 2, 3, 2, 1, 3, 2, 3, 2, 1],
          [2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2, 3, 2, 1]]
    l2 = [[2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2, 3],
          [2, 1, 2, 3, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 3],
          [1, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 3],
          [1, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2, 3],
          [2, 1, 2, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 3],
          [1, 2, 1, 3, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 3],
          [1, 3, 2, 1, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 2],
          [1, 3, 2, 3, 2, 1, 3, 2, 3, 2, 1, 3, 2, 3, 2]]
    return l1, l2

#################

def EquivalenceClasses(W):
    P = Poset((W,lambda p,q:p.bruhat_le(q)))
    E = [e for e in P.cover_relations() if not e[0].weak_le(e[1])]
    D = {tuple(e):[] for e in E}
    for w in CommutationClasses(W):
        L = WSortingOrder(W,w)
        for e in L.cover_relations():
            if not e[0].weak_le(e[1]):
                D[tuple(e)].append(tuple(w))
    return Stats(list(D),[lambda e:Set(D[e])],ret = True)

def Implications(W):
    P = Poset((W,lambda p,q:p.bruhat_le(q)))
    E = [e for e in P.cover_relations() if not e[0].weak_le(e[1])]
    D = {tuple(e):[] for e in E}
    for w in CommutationClasses(W):
        L = WSortingOrder(W,w)
        for e in L.cover_relations():
            if not e[0].weak_le(e[1]):
                D[tuple(e)].append(tuple(w))
    S = [s[0] for s in Stats(list(D),[lambda e:Set(D[e])],ret = True)]
    return Poset((range(len(S)),lambda i,j:S[i].issubset(S[j])))

def Dissectors_Test(W,w):
    L = WSortingOrder(W,w)
    for x in L.meet_irreducibles():
        P = L.subposet([y for y in L if not L.is_lequal(y,x)])
        m = P.minimal_elements()
        if len(m)!=1:
            print(x,m)

def JoinIrreducibles_Test(W,w):
    L = SortingOrder(W,w)
    P = L.join_irreducibles_poset()
    C,C2 = list(P)+[Set()],[]
    for i in range(len(C)):
        for j in range(i+1,len(C)):
            I = C[i].union(C[j])
            if I not in C+C2:
                C2.append(I)
    while C2!=[]:
        D = []
        for x in C:
            for y in C2:
                I = x.union(y)
                if all(I not in X for X in [C,C2,D]):
                    D.append(I)
        C.extend(C2)
        C2 = D
    return LatticePoset((C,lambda p,q:p.issubset(q))) == L

def AntiPart_to_JoinIrrPoset(P,L):
    """
    L = partition du poset P en antichaines
    """
    D = {k:i for i,l in enumerate(L) for k in l}
    return Poset(([Set([D[i] for i in P.order_ideal([k])]) for k in P], lambda p,q:p.issubset(q)))

def MeetIrreduciblesPoset(P,L):
    """
    L = partition du poset P en antichaines
    """
    return Poset((range(len(L)),lambda p,q:all(any(P.is_lequal(x,y) for x in L[p]) for y in L[q])))


def JoinClosure(L):
    C,C2 = list(L),[]
    for i in range(len(C)):
        for j in range(i+1,len(C)):
            I = C[i].union(C[j])
            if I not in C+C2:
                C2.append(I)
    while C2!=[]:
        D = []
        for x in C:
            for y in C2:
                I = x.union(y)
                if all(I not in X for X in [C,C2,D]):
                    D.append(I)
        C.extend(C2)
        C2 = D
    if Set() not in C:
        C.append(Set())
    return LatticePoset((C,lambda p,q:p.issubset(q)))


def AntiPartitions(P):
    if len(P)==0:
        return [[]]
    l = sorted(P)
    P2 = P.subposet([x for x in P if not P.is_lequal(x,l[0]) and not P.is_lequal(l[0],x)])
    LP = []
    for a in P2.antichains():
        a2 = [l[0]] + a
        P3 = P.subposet([x for x in P if x not in a2])
        for p in AntiPartitions(P3):
            LP.append([a2] + p)
    return LP

def ReducedAntiPartitions(P):
    return [p for p in AntiPartitions(P) if P.is_isomorphic(AntiPart_to_JoinIrrPoset(P,p))]

#énumération des partitions en antichaînes réduites
# à mettre en parallèle avec les partitions en antichaînes propres

def ChainProdPoset(l):
    """
    renvoie le poset [m]*[n]... pour l = [m,n]
    """
    return Poset((cartesian_product([range(k) for k in l]), lambda t1,t2 : all(i<=j for i,j in zip(t1,t2))))

def DisjointChainPoset(l):
    return Poset(([(k,i) for i,a in enumerate(l) for k in range(a)], lambda t1,t2 : t1[0]<=t2[0] and t1[1]==t2[1]))

def JoinDisLattices(n):
    """
    renvoie les treillis join-distributifs à n join-irréductibles
    """
    LL = []
    for P in Posets(n):
        Ll = []
        for p in AntiPartitions(P):
            J = AntiPart_to_JoinIrrPoset(P,p)
            if P.is_isomorphic(J):
                L = JoinClosure(J)
                if all(not M.is_isomorphic(L) for M in Ll):
                    Ll.append(L)
        LL.extend(Ll)
    return LL

"""
treillis à n join-irréductibles et k meet-irréductibles (k<=n)

n = 3 :
(3,) 5

n = 4
(3,) 1
(4,) 16

n = 5
(4,) 10
(5,) 63

n = 6
(4,) 7
(5,) 108
(6,) 318

n = 7
(4,) 1
(5,) 175
(6,) 1269
(7,) 2045
"""

############# Sorting orders maxflag ############

def MaxflagSortingOrders(W, tilings = False):
    if tilings:
        n = W.n
        if W.cartan_type()[0] != 'A':
            raise Exception('W must be of type A')
    l1, l2 = [], []
    for w in CommutationClasses(W):
        L = SortingOrder(W, w)
        if is_maxflag(L) == True:
            l1.append(w)
        else:
            l2.append(w)
    if tilings:
        return [W0ReducedWord(n, w) for w in l1], [W0ReducedWord(n, w) for w in l2]
    return l1, l2

def MinflagSortingOrders(W, tilings = False):
    if tilings:
        n = W.n
        if W.cartan_type()[0] != 'A':
            raise Exception('W must be of type A')
    l1, l2 = [], []
    for w in CommutationClasses(W):
        L = SortingOrder(W, w)
        if is_maxflag(L.dual()) == True:
            l1.append(w)
        else:
            l2.append(w)
    if tilings:
        return [W0ReducedWord(n, w) for w in l1], [W0ReducedWord(n, w) for w in l2]
    return l1, l2