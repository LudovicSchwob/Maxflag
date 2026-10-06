from Transformations_ASM import *

#Statistiques Alpha
def Alpha_ASM(M):
    if isinstance(M,ASM):
        M = M.main
    for k in range(len(M)):
        if M[k][0]==1:
            return k+1
def Alpha_Gog(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    for k in range(1,n):
        if G[-k-1][0]!=1:
            return k
    return n

#NE DONNE PAS LES BONNES VALEURS, POURQUOI ?
def Alpha_Magog2(M):
    S,n = Schutz_Tri(M).main,len(M.main)
    l = Alpha_Gog(S)
    if l!=n:
        for k in range(1,l):
            if S[n-l+k-1][k]==S[n-l+k][k]:
                return l-1
    return l

#Statistiques Beta
def Beta_ASM(M):
    if isinstance(M,ASM):
        M = M.main
    for k in range(len(M)):
        if M[-1][k]==1:
            return k+1
def Beta_Gog(G):
    if isinstance(G,GT_Tri):
        G = G.main
    return G[0][0]

def Beta_Magog(M):
    if isinstance(M,GT_Tri):
        M = M.main
    n,s = len(M),1
    for k in range(n-1):
        s += M[-1][k+1]-M[-2][k]
    return s
def Beta_Gogam(G):
    if isinstance(G,GT_Tri):
        G = G.main
    return G[0][0]

#Statistiques Gamma
def Gamma_ASM(M):
    if isinstance(M,ASM):
        M = M.main
    for k in range(len(M)):
        if M[k][-1]==1:
            return k+1
def Gamma_Gog(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    for k in range(1,n):
        if G[-k-1][-1]!=n:
            return k
    return n

def Gamma_Magog1(M):
    if isinstance(M,GT_Tri):
        M = M.main
    n,s = len(M),1
    for k in range(1,n):
        if M[k][-1] == k+1:
            s += 1
    return s

def Gamma_Magog(M):
    if isinstance(M,GT_Tri):
        M = M.main
    n= len(M)
    for i in range(n):
        if M[-i-1][-1]==n-i:
            k = n-i
            break
    i,j = k,k
    for l in range(1,k):
        if M[i-1][j-2]==M[i-2][j-2]:
            i,j = i,j-1
        else:
            i,j = i-1,j-1
    return i
def Gamma_Gogam(G):
    return Gamma_Magog(G.Schutz())

#Statistiques Delta
def Delta_ASM(M):
    if isinstance(M,ASM):
        M = M.main
    for k in range(len(M)):
        if M[0][k]==1:
            return k+1
def Delta_Gog(G):
    if isinstance(G,GT_Tri):
        G = G.main
    for k in range(1,len(G)):
        if G[-2][k-1]!=k:
            return k
    return len(G)

def Gamma_Magog2(M):
    B,x = M.to_Bool().main,0
    for k in range(len(B)):
        x += 1-B[-x-1][-k+x-1]
    return x


def Stat_Gog1(G):
    if isinstance(G,GT_Tri):
        G = G.main
    for k in range(len(G)):
        for i in range(k+1):
            if G[-k+i-1][i]!=i+1:
                return k
    return len(G)
def Stat_Magog1(M):
    if isinstance(M,GT_Tri):
        M = M.main
    for k in range(1,len(M)):
        if M[-1][k]!=1:
            return k
    return len(M)

def Stat_Gog2(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    def Stat_Aux(M):
        for k in range(1,len(M)):
            for i in range(k):
                if M[-k+i-1][-i-1]!=M[-k+i][-i-1]:
                    return k-1
        return len(M)
    for k in range(1,n):
        if Stat_Aux(Standard_Procedure([[G[-k+i-1][-i+j-1] for j in range(i+1)] for i in range(k+1)]))!=k+1:
            return k-1
    return n-1
def Stat_Magog2(M):
    n,M = len(M),GT_Tri(M).Schutz()
    if type(M)==GT_Tri:
        M = M.main
    for k in range(1,n):
        for i in range(k):
            if M[-k+i-1][-i-1]!=M[-k+i][-i-1]:
                return k-1
    return n-1

def Stat_Gog3(G):
    if type(G)==GT_Tri:
        G = G.main
    n = len(G)
    for k in range(1,n):
        for i in range(n-k):
            if G[-k][i+1]!=G[-k-1][i] and G[-k][i+1]!=G[-k][i]+1:
                return k-1
    return n-1
def Stat_Magog3(M):
    n,M = len(M),GT_Tri(M).Schutz()
    if type(M)==GT_Tri:
        M = M.main
    for k in range(1,n):
        for i in range(n-k):
            if M[-k][i+1]!=M[-k-1][i]:
                return k-1
    return n-1

def Gog_Dec_DS(G):
    if type(G)==GT_Tri:
        G = G.main
    l,n = [],len(G)
    for k in range(1,n):
        t = True
        for i in range(k):
            if G[-k+i-1][i]!=G[-k+i][i+1]:
                t = False
                break
        if t:
            l.append(k)
    return tuple(l)
def Magog_Dec_DS1(M):
    if type(M)==GT_Tri:
        M = M.main
    n,s = len(M),1
    for k in range(1,n):
        t = True
        for i in range(k):
            if M[-k+i-1][i]+1!=M[-k+i][i+1]:
                t = False
                break
        if t:
            s+=1
    return s

def Magog_Dec_DS2(M):
    if isinstance(M,GT_Tri):
        M = M.main
    l,n = [],len(M)
    for k in range(1,n):
        if M[-1][k]==k+1:
            if M[-1][:k]==M[k-1]:
                l.append(k)
    return tuple(l)

def Gog_Dec_SS(G):
    if isinstance(G,GT_Tri):
        G = G.main
    l,n = [],len(G)
    for k in range(1,n):
        t = True
        for i in range(n-k):
            if G[-k-1][i]!=i+1:
                t = False
                break
        if t:
            l.append(k)
    return tuple(l)
def Magog_Dec_SS2(M):
    if isinstance(M,GT_Tri):
        M = M.main
    l,n = [],len(M)
    for k in range(1,n):
        t = True
        for i in range(k):
            for j in range(min(i+1,n-k)):
                if M[k+j-1][i]!=M[k+j][i+1]:
                    t = False
                    break
            if not t:
                break
        if t:
            for i in range(k+1,n-1):
                for j in range(i-k):
                    if M[i][k+j+1]-M[i-j-1][k]>j+1:
                        t = False
                        break
                    if M[i][k+j+1]-M[i-j-1][k]<M[i+1][k+j+1]-M[i-j][k]:
                        t = False
                        break
                if not t:
                    break
            for j in range(n-1-k):
                if M[n-1][k+j+1]-M[n-j-2][k]>j+1:
                    t = False
                    break
            if t:
                l.append(k)
    return tuple(l)

def stat0(G):
    if isinstance(G,GT_Tri):
        G = G.main
    s = 0
    for k in range(len(G)):
        for i in range(k):
            if G[-k-1+i][i]==G[-k+i][i+1]:
                s+=1
    return s
def stat0bis(G):
    if isinstance(G,GT_Tri):
        G = G.main
    s = 0
    for k in range(len(G)):
        for i in range(k):
            if G[-k-1+i][i]==G[-k+i][i]:
                s+=1
    return s

def ASM_minus_1s(M):
    if isinstance(M,ASM):
        M = M.main
    n = len(M)-1
    return sum(1 for k in range(n) for i in range(1,len(M[0])) if M[k][i]==-1)

def Gog_minus_1s(G):
    if isinstance(G,GT_Tri):
        G = G.main
    return sum(1 for k in range(1,len(G)) for i in range(k) if G[k][i]!=G[k-1][i] and G[k][i+1]!=G[k-1][i])

def ASM_inv(M):
    if isinstance(M,ASM):
        M = M.main
    I,n,m,l,s = 0,len(M),len(M[0]),[],0
    if n==1:
        return 0
    for k in range(m-1):
        s += M[-1][k]
        I += s*M[-2][k+1]
        l.append(s)
    for k in range(2,n):
        s = 0
        for i in range(m-1):
            s += M[-k][i]
            l[i] += s
            I += M[-k-1][i+1]*l[i]
    return I
def ASM_coinv(M):
    if isinstance(M,ASM):
        M = M.main
    I,n,m,l,s = 0,len(M),len(M[0]),[],0
    if n==1:
        return 0
    for k in range(m-1):
        s += M[0][k]
        I += s*M[1][k+1]
        l.append(s)
    for k in range(2,n):
        s = 0
        for i in range(m-1):
            s += M[k-1][i]
            l[i] += s
            I += M[k][i+1]*l[i]
    return I
#!!!! ASM_inv <=> Gog_inv + # de -1 !!!!
def Gog_inv(G):
    if isinstance(G,GT_Tri):
        G = G.main
    s,n = 0,len(G)
    for k in range(1,n):
        for i in range(k):
            if G[k-1][i]==G[k][i]:
                s += 1
    return s
def Gog_coinv(G):
    if isinstance(G,GT_Tri):
        G = G.main
    s,n = 0,len(G)
    for k in range(1,n):
        for i in range(k):
            if G[k-1][i]==G[k][i+1]:
                s += 1
    return s

def DPP_parts(P):
    if type(P)==DPP:
        P = P.main
    return sum(len(P[k]) for k in range(len(P)))

def DPP_special_parts(P):
    """
    P doit être de type DPP
    """
    P = P.main
    return sum(1 for k in range(len(P)) for i in range(len(P[k])) if P[k][i]<=i)

def Delta_DPP(P):
    P,n = P.main,P.order
    if len(P)==0:
        return 1
    for k in range(len(P[0])):
        if P[0][k]!=n:
            return k+1
    return len(P[0])+1
def Beta_DPP(P):
    P,n = P.main,P.order
    s = int(len(P)>0 and len(P[0])==n-1)
    for k in range(min(len(P),2)):
        for i in P[k]:
            if i==n-1:
                s += 1
    return n-s

def Bool_des(B):
    B = B.main
    return sum(1 for b in B for k in range(1,len(b)) if b[k-1]>b[k])

def Gapless(T,d = 1,right = True):
    if isinstance(T,GT_Tri):
        T = T.main
    for k in range(1,len(T)):
        for i in range(k):
            if right:
                if T[k][i+1]-T[k-1][i]>d:
                    return False
            else:
                if T[k-1][i]-T[k][i]>d:
                    return False
    return True

def Gog_avoid132(G,side = "right"):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    if side == "left":
        for i in range(n-2):
            for j in range(i+1):
                if G[i][j]==G[i+1][j]>G[i+2][j]:
                    return False
    elif side == "right":
        for i in range(n-2):
            for j in range(i+1):
                if G[i][j]==G[i+1][j+1]<G[i+2][j+2]:
                    return False
    else:
        raise Exception("side must be 'left' or 'right'")
    return True

def BiGapless(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    for k in range(1,n):
        for i in range(k):
            if G[k][i+1]-G[k-1][i]>1 or G[k-1][i]-G[k][i]>1:
                return False
    return True
def RowGapless(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    for k in range(1,n):
        for i in range(k):
            if G[k][i+1]-G[k][i]>1:
                return False
    return True

def RowIncreasing(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    for k in range(1,n):
        for i in range(k):
            if G[k][i+1]<G[k][i]:
                return False
    return True
def DiagIncreasing(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n = len(G)
    for k in range(1,n):
        for i in range(k):
            if G[k][i+1]<G[k-1][i]:
                return False
    return True

def Non1Increasing(M):
    if isinstance(M,GT_Tri):
        M = M.main
    n = len(M)
    for k in range(1,n):
        for i in range(k):
            if M[k][i]==M[k][i+1]!=1:
                return False
    return True

def Stat_GogDownCover1(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n,l = len(G),[]
    for k in range(n-1):
        if G[k][0]>G[k+1][0]:
            l.append(k)
        for i in range(1,k+1):
            if G[k][i]>max(G[k-1][i-1],max(G[k+1][i],G[k][i-1]+1)):
                l.append(k-i)
    if l==[]:
        return -1
    return max(l)

def Stat_MagogDownCover1(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n,l = len(G),[]
    for k in range(n-1):
        for i in range(k):
            if G[k][i+1]>max(G[k-1][i],G[k+1][i+1]):
                l.append(k-i-1)
    for k in range(1,n):
        if G[-1][k]>G[-2][k-1]:
            l.append(n-k-1)
    if l==[]:
        return -1
    return max(l)

def Stat_GogDownCover2(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n,l = len(G),[]
    for k in range(n-1):
        if G[k][0]>G[k+1][0]:
            l.append(k+1)
        for i in range(1,k+1):
            if G[k][i]>max(G[k-1][i-1],max(G[k+1][i],G[k][i-1]+1)):
                l.append(k+1)
    if l==[]:
        return -1
    return max(l)

def Stat_MagogDownCover2(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n,l = len(G),[]
    for k in range(n-1):
        for i in range(k):
            if G[k][i+1]>max(G[k-1][i],G[k+1][i+1]):
                l.append(n-i-1)
    for k in range(1,n):
        if G[-1][k]>G[-2][k-1]:
            l.append(n-k)
    if l==[]:
        return -1
    return max(l)

def Stat_GogDownCover3(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n,l = len(G),[]
    for k in range(n-1):
        if G[k][0]>G[k+1][0]:
            l.append(G[k][0]+k)
        for i in range(1,k+1):
            if G[k][i]>max(G[k-1][i-1],max(G[k+1][i],G[k][i-1]+1)):
                l.append(G[k][i]+k-i)
    if l==[]:
        return -1
    return max(l)

def Stat_MagogDownCover3(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n,l = len(G),[]
    for k in range(n-1):
        for i in range(k):
            if G[k][i+1]>max(G[k-1][i],G[k+1][i+1]):
                l.append(n+G[k][i+1]-i-2)
    for k in range(1,n):
        if G[-1][k]>G[-2][k-1]:
            l.append(n+G[-1][k]-k-1)
    if l==[]:
        return -1
    return max(l)

def is_normal(T):
    T = T.main
    n,t = len(T),False
    if T[0][0]==0:
        t = True
    for k in range(1,n):
        s = False
        if T[k][0]!=0:
            s = True
        else:
            for i in range(k):
                if T[k-1][i]!=T[k][i+1]:
                    s = True
                    break
        if not s:
            t = True
        else:
            if t:
                return False
    return True

def U_Stat_Magog(M,k):
    if isinstance(M,GT_Tri):
        M = M.main
    n = len(M)
    return M[n-k][n-k]+sum(M[n-k][t]-M[n-k-1][t] for t in range(1,n-k))+k-1-sum(1 for t in range(n-k+1,n) if M[t][t]<t+1)

def U_Stat_Gog(G,k):
    if isinstance(G,GT_Tri):
        G = G.main
    if k==1:
        return G[0][0]
    return sum(G[k-1])-sum(G[k-2])+sum(1 for a in G[k-2] if a not in G[k-1])

#prend une ASM en argument et renvoie le chemin de dyck donné par ses excédances
def ExcNCP(M):
    if isinstance(M,ASM):
        M = M.main
    C,d = Corner_Sum(M),[0]
    for i in range(1,len(C)):
        d.append(2*C[-i-1][i-1]-1)
        d.append(2*C[-i-1][i])
    return tuple(d)

def Vexillary(G):
    if isinstance(G,GT_Tri):
        G = G.main
    n,E = len(G),[]
    for k in range(n-1):
        if G[k][0]>G[k+1][0]:
            E.append((k,G[k][0]))
        for i in range(1,k+1):
            if G[k][i]>max(G[k-1][i-1],max(G[k+1][i],G[k][i-1]+1)):
                E.append((k,G[k][i]))
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            if E[i][1]<E[j][1] and E[i][0]<E[j][0]:
                return False
            if E[i][1]>E[j][1] and E[i][0]>E[j][0]:
                return False
    return True
