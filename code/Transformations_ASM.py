from Classes import *

#Bijection entre triangles Gog et ASM
##def Gog_to_ASM(G):
##    n=len(G)
##    M=[n*[0] for _ in range(n)]
##    for i in range(n):
##        for j in range(n-i):
##            M[i][G[-i-1][j]-1]=1
##    for i in range(1,n):
##        for j in range(n):
##            M[i-1][j]-=M[i][j]
##    return M
##
##def ASM_to_Gog(M):
##    n = len(M)
##    l,G = n*[0],[]
##    for k in range(n):
##        g = []
##        for i in range(n):
##            l[i] += M[-k-1][i]
##            if l[i]:
##                g.append(i+1)
##        G.append(g)
##    return G
##
##def Rotate(P):
##    L=[]
##    if P!=[]:
##        for i in range(len(P[0])):
##            t=0
##            for p in P:
##                if len(p)<=i:
##                    break
##                t+=1
##            l=P[t-1][i]*[t]
##            for j in range(1,t):
##                l+=(P[t-j-1][i]-P[t-j][i])*[t-j]
##            L.append(l)
##    return L
##
###Triangle Magog -> TSSCPP
##def TSSCPP(M):
##    P,n=[],len(M)
##    for i in range(1,n):
##        p=[]
##        for j in range(1,n):
##            a=M[n-max(i,j)+min(i,j)-1][n-max(i,j)]
##            if a>1:
##                p.append(a-1)
##            else:
##                break
##        if p==[]:
##            break
##        else:
##            P.append(p)
##    T=[n*[2*n]+n*[n] for _ in range(n)]+[n*[n] for _ in range(n)]
##    for i in range(len(P)):
##        T[n+i]+=P[i]
##        for j in range(len(P[i])):
##            T[n-i-1][n-j-1]-=P[i][j]
##    P=Rotate(P)
##    for i in range(len(P)):
##        for j in range(len(P[i])):
##            T[n+i][j]+=P[i][j]
##            T[n-i-1][2*n-j-1]-=P[i][j]
##    P=Rotate(P)
##    for i in range(len(P)):
##        for j in range(len(P[i])):
##            T[i][n+j]+=P[i][j]
##            T[2*n-i-1][n-j-1]-=P[i][j]
##    return T
##
###Tableau Gog associé à une permutation
def Perm_to_Gog(p):
    return Gog([sorted(p[:k+1]) for k in range(len(p))])
def Gog_to_Perm(G):
    if isinstance(G,GT_Tri):
        G = G.main
    p,n = [G[0][0]],len(G)
    for i in range(1,n):
        s = False
        if G[i][0]!=G[i-1][0]:
            p.append(G[i][0])
            s = True
        for j in range(1,i):
            if G[i][j]!=G[i-1][j] and G[i][j]!=G[i-1][j-1]:
                if s:
                    raise Exception('Pas une permutation')
                p.append(G[i][j])
                s = True
        if G[i][-1]!=G[i-1][-1]:
            if s:
                raise Exception('Pas une permutation')
            p.append(G[i][-1])
    return p
        
def Perm_to_Magog(p):
    l = list(p)
    for k in range(len(p)):
        l[k] -= sum(1 for a in l[k+1:] if a<l[k])
    return Magog([sorted(l[:k+1]) for k in range(len(p))])
##
###Triangle G.T. -> Tableau de Young semi-standard
##def Tri_to_Tab(T):
##    if type(T)==GT_Tri:
##        T = T.main()
##    Y,n,z = [],len(T),0
##    while z<n and T[-1][z]==0:
##        z += 1
##    for k in range(n-z):
##        y = T[k][0]*[k+1]
##        for i in range(1,n-k):
##            y += (T[k+i][-k-1]-T[k+i-1][-k-1])*[k+i+1]
##        Y.append(y)
##    return Y
##
##def Tab_to_Tri(Y):
##    T = []
##    for k in range(max(y[-1] for y in Y)):
##        t = []
##        for i in range(k+1):
##            if k-i<len(Y):
##                for j in range(len(Y[k-i])):
##                    if Y[k-i][j]>k+1:
##                        j -= 1
##                        break
##                t.append(j+1)
##            else:
##                t.append(0)
##        T.append(t)
##    return T
##                
##
###involution de Schützenberger sur les tableaux de Young
##def Schutz_Tab(T):
##    if type(T)==Young_Tab:
##        T = T.main()
##    T2 = [t[:] for t in T]
##    n = sum(len(t) for t in T)
##    m = max([t[-1] for t in T])+1
##    for k in range(n):
##        x,y=0,0
##        v = T2[0][0]
##        while True:
##            if y<len(T[x])-1 and T2[x][y+1]>0:
##                if x<len(T)-1 and y<len(T[x+1]) and T2[x+1][y]>0 and T2[x][y+1]>=T2[x+1][y]:
##                    T2[x][y] = T2[x+1][y]
##                    x+=1
##                else:
##                    T2[x][y] = T2[x][y+1]
##                    y+=1
##            else:
##                if x<len(T)-1 and y<len(T[x+1]) and T2[x+1][y]>0:
##                    T2[x][y] = T2[x+1][y]
##                    x+=1
##                else:
##                    T2[x][y] = -v
##                    break
##    for i in range(len(T)):
##        for j in range(len(T[i])):
##            T2[i][j]+=m
##    return Young_Tab(T2)
##
##def omega(T,j):
##    T[0][0] = T[1][0]+T[1][1]-T[0][0]
##    for k in range(1,j):
##        T[k][0] = T[k+1][0]+min(T[k+1][1],T[k-1][0])-T[k][0]
##        for i in range(1,k):
##            T[k][i] = max(T[k+1][i],T[k-1][i-1])+min(T[k+1][i+1],T[k-1][i])-T[k][i]
##        T[k][-1] = max(T[k+1][-2],T[k-1][-1])+T[k+1][-1]-T[k][-1]
##    return T
##
###involution de Schützenberger sur les triangles de Gelfand-Tsetlin
##def Schutz_Tri(T):
##    if type(T)==GT_Tri:
##        T = T.main()
##    n = len(T)
##    for k in range(1,n):
##        T = omega(T,n-k)
##    return GT_Tri(T)

def Gog_Direct_Sum(G1,G2):
    G1,G2 = G1.main,G2.main
    n1,n2 = len(G1),len(G2)
    return GT_Tri([[i+n1 for i in g] for g in G2]+[g+list(range(n1+1,n1+n2+1)) for g in G1])

def Gog_Skew_Sum(G1,G2):
    G1,G2 = G1.main,G2.main
    n1,n2 = len(G1),len(G2)
    return GT_Tri(G2+[list(range(1,n2+1))+[i+n2 for i in g] for g in G1])

##def Magog_Direct_Sum1(M1,M2):
##    M1,M2 = M1.main(),M2.main()
##    n1,n2 = len(M1),len(M2)
##    S = []
##    for k in range(n1+n2):
##        l = []
##        if k>=n2:
##            l += M1[k-n2]
##        l+=list(range(max(1,k-n2+2),min(n1+1,k+2)))
##        if k>=n1:
##            l += [i+n1 for i in M2[k-n1]]
##        S.append(l)
##    return GT_Tri(S)

def GT_Transpose(M):
    if type(M)==GT_Tri:
        M = M.main
    return [[M[-i+j][j] for j in range(i)] for i in range(1,len(M)+1)]

##def Magog_Skew_Sum1(M1,M2):
##    M1,M2 = M1.main(),M2.main()
##    S = []
##    n1,n2 = len(M1),len(M2)
##    M1,M2 = n2*[[]]+GT_Transpose(M1),n1*[[]]+GT_Transpose(M2)
##    for k in range(n1+n2):
##        while M1[k]!=[] and M1[k][0]==1:
##            del M1[k][0]
##        while M2[k]!=[] and M2[k][0]==1:
##            del M2[k][0]  
##    for k in range(n1+n2):
##        S.append((k+1-len(M2[k])-len(M1[k]))*[1]+M2[k]+[M1[k][i]+n2 for i in range(len(M1[k]))])
##    return GT_Tri(GT_Transpose(S))

def Magog_Direct_Sum(M1,M2):
    M1,M2 = M1.main,M2.main
    n1,n2 = len(M1),len(M2)
    return GT_Tri(M1+[M1[-1]+[i+n1 for i in g] for g in M2])

def Magog_Skew_Sum(M1,M2):
    M1,M2 = M1.main,M2.main
    return GT_Tri(M1+[sorted(M1[-1]+m) for m in M2])

def Gogam_Direct_Sum(G1,G2):
    G1,G2 = G1.main,G2.main
    n1,n2 = len(G1),len(G2)
    return GT_Tri([[i+n1 for i in g] for g in G2]+[g+[i+n1 for i in G2[-1]] for g in G1])
def Gogam_Skew_Sum(G1,G2):
    G1,G2 = G1.main,G2.main
    return GT_Tri(G2+[sorted(G2[-1]+m) for m in G1])

def Bool_Direct_Sum(B1,B2):
    B1,B2 = B1.main,B2.main
    n1 = len(B1)
    return Bool(B1+[(n1+1)*[0]]+[(n1+1)*[0]+b for b in B2])
def Bool_Skew_Sum(B1,B2):
    B1,B2 = B1.main,B2.main
    n1 = len(B1)
    return Bool(B1+[(n1+1)*[1]]+[b+(n1+1)*[1] for b in B2])

##def Magog_Skew_Sum2(M1,M2):
##    M1,M2 = M1.main(),M2.main()
##    n1,n2 = len(M1),len(M2)
##    return GT_Tri(M1+[min(k+1,n1)*[1]+M1[-1][:max(0,n1-k-1)]+[M2[k][i]+M1[-1][max(n1+i-k-1,0)]-1 for i in range(k+1)] for k in range(n2)])

##def Magog_Direct_Sum3(M1,M2):
##    M1,M2 = M1.main(),M2.main()
##    n1,n2 = len(M1),len(M2)
##    S = []
##    for k in range(n1+n2):
##        l = []
##        if k>=n2:
##            l += M1[k-n2]
##        l+=[M1[i][-1] for i in range(max(0,k-n2+1),min(n1,k+1))]
##        if k>=n1:
##            l += [i+n1 for i in M2[k-n1]]
##        S.append(l)
##    return GT_Tri(S)

##def Magog_Skew_Sum3(M1,M2):
##    M1,M2 = M1.main(),M2.main()
##    S = []
##    n1,n2 = len(M1),len(M2)
##    M1,M2 = n2*[[]]+GT_Transpose(M1),n1*[[]]+GT_Transpose(M2)
##    for k in range(n1+n2):
##        while M1[k]!=[] and M1[k][0]==1:
##            del M1[k][0]
##    for k in range(n1+n2):
##        if len(M2[k])>0:
##            a = M2[k][-1]
##        else:
##            a = 1
##        S.append(M2[k]+(k+1-len(M1[k])-len(M2[k]))*[a]+[M1[k][i]+n2 for i in range(len(M1[k]))])
##    return GT_Tri(GT_Transpose(S))

###chemin surdiagonal -> sous-diagonal
##def Sup_to_Sub(d):
##    n=len(d)
##    if n<=1:
##        return d
##    a=0
##    for k in range(1,n):
##        if d[-k-1]==n-k:
##            a=n-k
##            break
##    l=d[a:-1]
##    for k in range(n-a-1):
##        l[k]-=a+1
##    l,m=Sup_to_Sub(l),Sup_to_Sub(d[:a])
##    for k in range(n-a-1):
##        l[k]+=a+1
##    return [1]+m+l
###chemin sous-diagonal -> sur-diagonal
##def Sub_to_Sup(d):
##    n=len(d)
##    if n<=1:
##        return d
##    a=n
##    for k in range(1,n):
##        if d[k]==k+1:
##            a=k
##            break
##    l=d[a:]
##    for k in range(n-a):
##        l[k]-=a
##    l,m=Sub_to_Sup(l),Sub_to_Sup(d[1:a])
##    for k in range(n-a):
##        l[k]+=a
##    return m+l+[n]
##
##
##def Standard_Procedure(G):
##    n,G = len(G),[list(g[:]) for g in G]
##    for k in range(1,n):
##        for i in range(k):
##            if G[-i-1][-k-1+i]==G[-i-2][-k+i]:
##                for j in range(i+1):
##                    G[-i-1+j][-k+i] -= 1
##    return G
##
##def Reverse_Standard_Procedure(G):
##    n,G = len(G),[g[:] for g in G]
##    for k in range(n-1):
##        for i in range(1,n-k):
##            if G[k+i-1][k]==G[k+i][k]:
##                for j in range(n-k-i):
##                    G[k+i+j][k+j+1] += 1
##    return G

def Sym_DPP(P):
    N,P = P.order,P.main
    n = len(P)
    shape = list(range(1,N-n)) + [N-n+i-len(P[-i-1]) for i in range(n)]
    values = [[] for _ in range(N-1)]
    sp_values = [[] for _ in range(N-1)]
    for k in range(n):
        for i in range(len(P[k])):
            if P[k][i]<=i:
                values[i-1].append(N-i-k-2)
                sp_values[-k-1].append(i+1-P[k][i])
            else:
                values[P[k][i]-2].append(N-i-k-2)
    remp = [[] for _ in range(N-1)]
    for k in range(1,N):
        c = 0
        for a in range(N-1):
            if a not in values[-k] and len(remp[a])<shape[a]:
                remp[a].append(N-k+1)
                c += 1
            if c==N-k-len(values[-k]):
                break
    for k in range(N-1):
        remp[k] += [sp_values[k][-i-1] for i in range(len(sp_values[k]))]
    sym = [[] for _ in range(max(len(r) for r in remp))]
    for r in remp:
        for i in range(len(r)):
            sym[i].append(r[i])
    return DPP(sym,N)

def Corner_Sum(M):
    if type(M)!=list:
        M = M.main
    n, m = len(M),len(M[0])
    S = [l+[0] for l in M]+[(m+1)*[0]]
    for i in range(1,n):
        for j in range(m):
            S[-i-2][j] += S[-i-1][j]
    for i in range(n):
        for j in range(1,m):
            S[i][-j-2] += S[i][-j-1]
    return S

def Height_Function(M):
    if type(M)==ASM:
        M = M.main
    n = len(M)
    S = [m+[0] for m in M]+[(n+1)*[0]]
    for i in range(1,n):
        for j in range(n):
            S[-i-2][j] += S[-i-1][j]
    for i in range(n):
        for j in range(1,n):
            S[i][-j-2] += S[i][-j-1]
    return [[2*S[i][j]+i+j-n for j in range(n+1)] for i in range(n+1)]

def omega_Magog(T,j):
    n = len(T)
    for k in range(1,j):
        for i in range(1,k):
            T[k][i] = max(T[k+1][i],T[k-1][i-1])+min(T[k+1][i+1],T[k-1][i])-T[k][i]
        T[k][-1] = max(T[k+1][-2],T[k-1][-1])+min(T[k+1][-1],k+1)-T[k][-1]
    return T

#involution de Schützenberger sur les triangles de Gelfand-Tsetlin
def Schutz_Magog(T):
    if isinstance(T,GT_Tri):
        T = T.main
    T = [t[:] for t in T]
    n = len(T)
    for k in range(1,n):
        T = omega_Magog(T,n-k)
    return GT_Tri(T)

def omega_Gog(T,j):
    n = len(T)
    for k in range(1,j):
        T[k][0] = T[k+1][0] + min([T[k+1][1],T[k-1][0],T[k][1]-1]) - T[k][0]
        for i in range(1,k):
            T[k][i] = max([T[k+1][i],T[k-1][i-1],T[k][i-1]+1])+min([T[k+1][i+1],T[k-1][i],T[k][i+1]-1])-T[k][i]
        T[k][-1] = max([T[k+1][-2],T[k-1][-1],T[k][-2]])+T[k+1][-1]-T[k][-1]
    return T
def Schutz_Gog(T):
    if isinstance(T,GT_Tri):
        T = T.main
    T = [t[:] for t in T]
    n = len(T)
    for k in range(1,n):
        T = omega_Gog(T,n-k)
    return GT_Tri(T)
