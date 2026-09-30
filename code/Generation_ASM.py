from Transformations_ASM import *

def Gogs(n):
    L=[[list(range(1,n+1))]]
    for i in range(1,n):
        M = []
        for l in L:
            for k in range(l[0][0],l[0][1]+1):
                M.append([[k]]+l)
        L = M
        for j in range(1,n-i):
            M=[]
            for l in L:
                for k in range(max(l[1][j],l[0][j-1]+1),l[1][j+1]+1):
                    M.append([l[0]+[k]]+l[1:])
            L = M
    return [Gog(G) for G in L]

def Magogs(n):
    L=[[[1] for _ in range(n)]]
    for i in range(1,n):
        M = []
        for l in L:
            for k in range(l[-2][i-1],i+2):
                M.append(l[:-1] + [l[-1]+[k]])
        L=M
        for j in range(1,n-i):
            M=[]
            for l in L:
                for k in range(max(l[-2-j][i-1],l[-j][i]),i+2):
                    M.append(l[:-j-1]+[l[-j-1]+[k]]+l[-j:])
            L = M
    return [Magog(M) for M in L]

def CatalanTriangles(n):
    L=[[list(range(1,n+1))]]
    for i in range(1,n):
        M = []
        for l in L:
            for k in range(max(l[0][0],l[0][1]-1),min(l[0][0]+1,l[0][1])+1):
                M.append([[k]]+l)
        L = M
        for j in range(1,n-i):
            M=[]
            for l in L:
                for k in range(max(l[1][j],l[1][j+1]-1),min(l[1][j]+1,l[1][j+1])+1):
                    M.append([l[0]+[k]]+l[1:])
            L = M
    return [CatalanTriangle(G) for G in L]

def TSSCPPs(n):
    return [M.to_TSSCPP() for M in Magogs(n)]
def ASMs(n):
    return [G.to_ASM() for G in Gogs(n)]

def Gogams(n):
    return [M.Schutz() for M in Magogs(n)]

def Irreduc_Gogs(n):
    I,I1,I2,M = [[GT_Tri([[1]])]],[[GT_Tri([[1]])]],[[GT_Tri([[1]])]],[[GT_Tri([[1]])]]
    for k in range(1,n):
        II = Gogs(k+1)
        II1,II2 = list(II),list(II)
        M.append(list(II))
        for i in range(k):
            for ii in I1[i]:
                for jj in M[k-i-1]:
                    ds = Gog_Direct_Sum(ii,jj)
                    del II1[II1.index(ds)]
                    del II[II.index(ds)]
            for ii in I2[i]:
                for jj in M[k-i-1]:
                    ds = Gog_Skew_Sum(ii,jj)
                    del II2[II2.index(ds)]
                    del II[II.index(ds)]
        I.append(II)
        I1.append(II1)
        I2.append(II2)
    return I[-1]

def Irreducs(n,Objects,Direct_Sum,Skew_Sum):
    I,I1,I2,M = [Objects(1)],[Objects(1)],[Objects(1)],[Objects(1)]
    for k in range(1,n):
        II = Objects(k+1)
        II1,II2 = list(II),list(II)
        M.append(list(II))
        for i in range(k):
            for ii in I1[i]:
                for jj in M[k-i-1]:
                    ds = Direct_Sum(ii,jj)
                    del II1[II1.index(ds)]
                    del II[II.index(ds)]
            for ii in I2[i]:
                for jj in M[k-i-1]:
                    ds = Skew_Sum(ii,jj)
                    del II2[II2.index(ds)]
                    del II[II.index(ds)]
        I.append(II)
        I1.append(II1)
        I2.append(II2)
    return I[-1]

def Irreduc_Magogs(n):
    return Irreducs(n,Magogs,Magog_Direct_Sum,Magog_Skew_Sum)
def Irreduc_Gogams(n):
    return Irreducs(n,Gogams,Gogam_Direct_Sum,Gogam_Skew_Sum)
def Irreduc_Bools(n):
    return Irreducs(n,Bools,Bool_Direct_Sum,Bool_Skew_Sum)

def Gog_Pentagons(n,k,l,m):
    if m>n:
        m = n
    if l>m:
        l = m
    if k>m:
        k = m
    P = []
    for i in range(m-l,n-l+1):
        P.append([(m-l)*[None]+[i+1]])
    for i in range(1,k+l-m):
        Q = []
        for p in P:
            for j in range(p[0][-1],n-l+i+1):
                q = [r[:] for r in p]
                q[0].append(j+1)
                Q.append(q)
        P = Q
    P = [[r[:]+(m-k)*[None] for r in p] for p in P]
    for i in range(1,m):
        Q = []
        for p in P:
            if m-l-i>0:
                Q.append([[None]]+p)
            else:
                if m-l-i==0:
                    s = 1
                else:
                    s = p[0][0]
                for a in range(s,p[0][1]+1):
                    Q.append([[a]]+p)
        P = Q
        for j in range(1,m-i):
            Q = []
            for p in P:
                if m-l-i-j>0 or j-k+1>0:
                    Q.append([p[0]+[None]]+p[1:])
                else:
                    if m-l-i-j==0:
                        s = j+1
                    else:
                        s = max(p[1][j],p[0][j-1]+1)
                    if j-k+1==0:
                        t = n-m+i+j+1
                    else:
                        t = p[1][j+1]
                    for a in range(s,t+1):
                        Q.append([p[0]+[a]]+p[1:])
            P = Q
    return [GT_Tri(G) for G in P]

def DPPs(n):
    if n==1:
        return [DPP([],1)]
    M = [[[k]] for k in range(2,n+1)]
    L = [[]] + M
    while M!=[]:
        M2 = []
        for l in M:
            nl = len(l[-1])
            if l[-1][0]<n:
                if len(l)==1:
                    for i in range(l[0][0]-nl):
                        M2.append([[k+1 for k in l[0]]+(i+1)*[1]])
                elif nl<len(l[-2])-1 and l[-1][0]<len(l[-2]):
                    s,m = True,[]
                    for k in range(nl):
                        if l[-1][k]+1==l[-2][k+1]:
                            s = False
                            break
                        m.append(l[-1][k]+1)
                    if s:
                        k = nl+1
                        while k<len(l[-2]) and l[-2][k]>1:
                            k += 1
                        for i in range(1,min(k,l[-1][0]+1)-nl):
                            M2.append(l[:-1]+[m+i*[1]])
            if len(l[-1])>1:
                if len(l)==1:
                    if l[0][-1]<l[0][-2]:
                        M2.append([l[0][:-1]+[l[0][-1]+1]])
                elif l[-1][-1]<l[-1][-2] and l[-1][-1]<l[-2][nl]-1:
                    M2.append(l[:-1]+[l[-1][:-1]+[l[-1][-1]+1]])
                for k in range(2,min(nl+1,l[-1][1])):
                    M2.append(l+[[k]])
        M = M2
        L += M
    return [DPP(P,n) for P in L]

def Bools(n):
    if n==1:
        return [Bool([])]
    L = [[[0], [0]], [[1], [1]]]
    for k in range(2,n):
        M = []
        for l in L:
            M.append([l[0]+[0]]+l[1:]+[[0]])
            M.append([[l[0][0]+1]+l[0][1:]+[0]]+l[1:]+[[1]])
        L = M
        for i in range(1,k):
            M = []
            for l in L:
                ll = l[0][:]
                ll[i] += 1
                M.append([ll]+l[1:-1]+[l[-1]+[1]])
                if l[0][i-1]<=l[0][i]+1:
                    M.append(l[:-1]+[l[-1]+[0]])
            L = M
    return [Bool(l[1:]) for l in L]

def ASTs(n):
    L = [[[k]] for k in range(1,2*n)]
    for k in range(1,n):
        M = []
        for l in L:
            b = l[-1][0]
            for i in range(k+1,min(b+1,2*n-k)):
                M.append(l+[[i]])
            if b<=k:
                M.append(l+[[b]])
        L = M
        for i in range(1,k):
            M = []
            for l in L:
                a,b = l[-2][i-1],l[-2][i]
                if a==l[-1][-1]:
                    for j in range(max(k,a)+1,min(b+1,2*n-k)):
                        M.append(l[:-1]+[l[-1]+[j]])
                else:
                    for j in range(max(k+1,a),min(b+1,2*n-k)):
                        M.append(l[:-1]+[l[-1]+[j]])
                    if a>=2*n-k:
                        M.append(l[:-1]+[l[-1]+[a]])
                if b<=k:
                    M.append(l[:-1]+[l[-1]+[b]])
            L = M
        M = []
        for l in L:
            a = l[-2][-1]
            if a==l[-1][-1]:
                for j in range(max(k,a)+1,2*n-k):
                    M.append(l[:-1]+[l[-1]+[j]])
            else:
                for j in range(max(k+1,a),2*n-k):
                    M.append(l[:-1]+[l[-1]+[j]])
                if a>=2*n-k:
                    M.append(l[:-1]+[l[-1]+[a]])
        L = M
    return [AST(T) for T in L]

def GT_Tris(l):
    L,n=[[l]],len(l)
    for i in range(1,n):
        M = []
        for l in L:
            for k in range(l[0][0],l[0][1]+1):
                M.append([[k]]+l)
        L = M
        for j in range(1,n-i):
            M=[]
            for l in L:
                for k in range(max(l[1][j],l[0][j-1]),l[1][j+1]+1):
                    M.append([l[0]+[k]]+l[1:])
            L = M
    return [GT_Tri(l) for l in L]


def Subpaths(n):
    L = [[1]]
    for k in range(1,n):
        M = []
        for l in L:
            for i in range(l[-1],k+2):
                M.append(l+[i])
        L = M
    return [Subpath(l) for l in L]
def Parkings(n):
    L = [[1]]
    for k in range(1,n):
        M = []
        for l in L:
            for i in range(1,k+2):
                M.append(l+[i])
        L = M
    return [Parking(l) for l in L]
def Paths(n,m):
    L = [[k] for k in range(0,m+1)]
    for _ in range(1,n):
        M = []
        for l in L:
            for i in range(l[-1],m+1):
                M.append(l+[i])
        L = M
    return [Path(l,m) for l in L]

def GenMagogs1(n):
    L = [n*[[1]]]
    for k in range(1,n):
        P = [[i] for i in range(1,k+2)]
        for _ in range(1,n-k):
            M = []
            for l in P:
                for i in range(1,l[-1]+1):
                    M.append(l+[i])
            P = M
        M = []
        for l in P:
            for m in L:
                m2 = [a.copy() for a in m]
                for i,a in enumerate(l):
                    m2[i+k].append(a)
                M.append(m2)
        L = M
    return [GT_Tri(l) for l in L]

def GenMagogs2(n):
    L = [[[1]]]
    for k in range(1,n):
        M = []
        for l in L[-1]:
            for i in range(l[-1],k+2):
                M.append(l+[i])
        L.append(M)
    L2 = [[[1]]]
    for k in range(1,n):
        M = []
        for l in L[k]:
            for m in L2:
                M.append(m+[l])
        L2 = M
    return [GT_Tri(l) for l in L2]
