def ASM_quotient(M):
    G = Gog(M.to_Tri())
    n = len(G.main)
    arc_list = []
    for m in G.canonical_meetands(form = 'irr'):
        a,b,c = m.main
        for s in Subsets(range(a+2,c), b-a-1):
            s2 = [i for i in range(a+2, c) if i not in s]
            arc_list.append((list(s), a+1, c, s2))
    return WOQuotient.fromArcList(arc_list, n)

def Maxflag_ASMs(n):
    l = []
    for M in ASMs(n):
        Q = ASM_quotient(M)
        L = LatticePoset(Q.poset())
#        print(Q)
#        print(Integer(len(L)), Integer(len(L)).factor())
        if is_maxflag(L) == True:
            l.append(M)
    return l

def Maxflag_Permutations(n):
    l = []
    for p in Permutations(n):
        M = Perm_to_Gog(p).to_ASM()
        Q = ASM_quotient(M)
        L = LatticePoset(Q.poset())
#        print(Q)
#        print(Integer(len(L)), Integer(len(L)).factor())
        if is_maxflag(L) == True:
            l.append(p)
    return l