from Generation_ASM import *
from Statistiques_ASM import *
import random
import string

def Draw_List(L,title = None,dim = (10,5),size = 70,style = None,save = False,show = True):
    x,y = dim
    if save and title==None:
        title,letters = '',string.ascii_letters
        for _ in range(10):
            title += random.choice(letters)
        print(f'File saved as {title}_.svg')
    for k in range((len(L)-1)//(x*y)+1):
        fig, ax = plt.subplots(y,x,figsize=(16,8),subplot_kw={'aspect': 'equal'})
        for i in range(y):
            for j in range(x):
                ax[i, j].axis('off')
                if k*x*y+i*x+j<len(L) and L[k*x*y+i*x+j]!=None:
                    L[k*x*y+i*x+j].draw(ax[i, j],size,style)
        if save:
            plt.savefig(title+str(k)+'.svg')
        if show:
            fig.suptitle(title)
            plt.show()

def Draw_Bij(L,bijs = [],title = '',dim = (10,2),size = 70,func = lambda x:x):
    if type(L)==dict:
        D,L = L,list(L)
        bijs = [lambda t:D[t]]
    b = len(bijs)
    x,y = dim
    for k in range((len(L)-1)//(x*y)+1):
        fig, ax = plt.subplots(y*(b+2)-1,x, figsize=(16,8),subplot_kw={'aspect': 'equal'})
        for j in range(x):
            for i in range(y*(b+2)-1):
                ax[i, j].axis('off')
            for i in range(y):
                if k*x*y+y*j+i<len(L) and L[k*x*y+y*j+i]!=None:
                    func(L[k*x*y+y*j+i]).draw(ax[(b+2)*i, j],size)
                    for z in range(b):
                        T = bijs[z](L[k*x*y+y*j+i])
                        if T!=None:
                            func(T).draw(ax[(b+2)*i+z+1, j],size)
        fig.suptitle(title)
        plt.show()

def Stats(T,stats,ret = False):
    D = {}
    for t in T:
        s = []
        for stat in stats:
            s.append(stat(t))
        s = tuple(s)
        if s in D:
            D[s].append(t)
        else:
            D[s] = [t]
    if ret:
        return D
    else:
        for s in sorted(list(D)):
            print(s,len(D[s]))

def Draw_Stat(T,stats,bijs = [],func = lambda x:x,style = None,one = False):
    D = {}
    for t in T:
        s = []
        for stat in stats:
            s.append(stat(t))
        s = tuple(s)
        if s in D:
            D[s].append(func(t))
        else:
            D[s] = [func(t)]
    if one:
        for s in list(D):
            if len(D[s])>1:
                del D[s]
    statname = str(tuple(stat.__name__ for stat in stats))
    if bijs:
        for s in D:
            Draw_Bij(D[s],bijs,statname+' = '+str(s))
    else:
        for s in D:
            Draw_List(D[s],statname+' = '+str(s),style = style)
def prod_aux(m,n):
    if m==0:
        return [[]],[list(range(n))]
    if n==0:
        return [list(range(m))],[[]]
    P,Q,a = prod_aux(m-1,n),prod_aux(m,n-1),n+m-1
    return Q[0]+[p+[a] for p in P[0]],[q+[a] for q in Q[1]]+P[1]

def ASM_Product(A,B):
    A,B = A.main(),B.main()
    n,m,L = len(A),len(B),[]
    P,Q = prod_aux(n,m)
    for k in range(len(P)):
        M = [(m+n)*[0] for _ in range(n+m)]
        for i in range(n):
            for j in range(n):
                M[j][P[k][i]] = A[j][i]
        for i in range(m):
            for j in range(m):
                M[j+n][Q[k][i]] = B[j][i]
        L.append(ASM(M))
    return L

def Gog_Product(G1,G2):
    return [A.to_Tri() for A in ASM_Product(G1.to_ASM(),G2.to_ASM())]

#G = gog associé à une permutation
def Gog_to_Magog(G):
    if isinstance(G,GT_Tri):
        G = G.main
    l,n = [G[0][0]],len(G)
    for k in range(1,n):
        for i in range(k+1):
            if (i==0 or G[k][i]!=G[k-1][i-1]) and (i==k or G[k][i]!=G[k-1][i]):
                l.append(G[k][i])
                break
    for k in range(1,n):
        l[-k] -= sum(1 for a in l[:-k] if l[-k]>a)
    return Magog([sorted(l[k:]) for k in range(n-1,-1,-1)])
def Magog_to_Gog(G):
    if isinstance(G,GT_Tri):
        G = G.main
    l,n = [1],len(G)
    for k in range(1,n):
        for i in range(k,-1,-1):
            if i==0 or G[k][i]!=G[k-1][i-1]:
                l.append(G[k][i])
                break
    p = n*[0]
    for k in range(n):
        t = l[-k-1]
        t += sum(p[:t])
        while p[t-1]:
            t += 1
        p[t-1] = 1
        l[-k-1] = t
    return Gog([sorted(l[k:]) for k in range(n-1,-1,-1)])

#G triangle de GT avec 1ère ligne 123...n
def GT_to_Magog(G):
    if isinstance(G,GT_Tri):
        G = G.main
    M,n = [],len(G)
    for k in range(n):
        m = (n+1-G[-k-1][-1])*[1]
        for i in range(k):
            m += (G[-k+i-1][-i-1]-G[-k+i][-i-2])*[i+2]
        M.append(m)
    return Magog(M)
def Magog_to_GT(M):
    if isinstance(M,GT_Tri):
        M = M.main
    n = len(M)
    G = [[] for _ in range(n)]
    for k in range(n):
        m = (n+1-k-M[-k-1][-1])*[k+1]
        for i in range(1,n-k):
            m += (M[-k-1][-i]-M[-k-1][-i-1])*[k+i+1]
        for i in range(n-k):
            G[-i-1].append(m[i])
    return GT_Tri(G)

#a,b = ± 1
# symétrie d'axe horizontal : -1,1,False
# symétrie d'axe vertical : 1,-1,False
#transposée : 1,1,True
def Sym_ASM(M,a = 1,b = 1,ref = False):
    M = M.main
    n,m = len(M),len(M[0])
    d1,d2 = (1-a)//2,(1-b)//2
    M = [[M[a*i-d1][b*j-d2] for j in range(m)] for i in range(n)]
    if ref:
        M = [[M[j][i] for j in range(n)] for i in range(m)]
    return ASM(M)

def Cycle_type(B):
    L,C = list(B),[]
    while L!=[]:
        s,l0 = 1,L[0]
        del L[0]
        l = B[l0]
        while l!=l0:
            del L[L.index(l)]
            l = B[l]
            s += 1
        C.append(s)
    return sorted(C)

#p et q chemins sous-diagonaux
def Triangles(p,q,magog = False):
    n = len(p)
    L = [[[k] for k in reversed(q[1:])]+[p]]
    for k in range(2,n):
        for i in range(n-k):
            M = []
            if magog:
                b = i+3
            else:
                b = k+i+2
            for l in L:
                for j in range(max(l[-k][i],l[-k+1][i+1]),b):
                    m = [a.copy() for a in l]
                    m[-k].append(j)
                    M.append(m)
            L = M
    return [GT_Tri(l) for l in L]
