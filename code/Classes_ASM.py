import matplotlib.pyplot as plt
from math import *
from Transformations_ASM import *
import random
import string

tfont = {'fontname':'times'}

def Rotate(P):
    L=[]
    if P!=[]:
        for i in range(len(P[0])):
            t=0
            for p in P:
                if len(p)<=i:
                    break
                t+=1
            l=P[t-1][i]*[t]
            for j in range(1,t):
                l+=(P[t-j-1][i]-P[t-j][i])*[t-j]
            L.append(l)
    return L
def omega(T,j):
    T[0][0] = T[1][0]+T[1][1]-T[0][0]
    for k in range(1,j):
        T[k][0] = T[k+1][0]+min(T[k+1][1],T[k-1][0])-T[k][0]
        for i in range(1,k):
            T[k][i] = max(T[k+1][i],T[k-1][i-1])+min(T[k+1][i+1],T[k-1][i])-T[k][i]
        T[k][-1] = max(T[k+1][-2],T[k-1][-1])+T[k+1][-1]-T[k][-1]
    return T
def omega_ASM(M,j):
    n,m = len(M),len(M[0])
    s = m*[0]
    for k in range(j):
        l1,l2 = list(M[-k-2]),list(M[-k-1])
        for i in range(m):
            if s[i]==0 and l1[i]==-1:
                t1 = i-1
                while M[-k-2][t1]==0 and s[t1]==0:
                    t1 -= 1
                l1[t1] -= 1
                l2[t1] += 1
                t2 = i+1
                while M[-k-2][t2]==0 and s[t2]==0:
                    t2 += 1
                l1[t2] -= 1
                l2[t2] += 1
                l1[i],l2[i] = 0,0
                l1[t1+t2-i],l2[t1+t2-i] = 1,-1
        M[-k-1],M[-k-2] = l1,l2
        for i in range(m):
            s[i] += l1[i]
    return M
def GT_Transpose(M):
    return [[M[-i+j][j] for j in range(i)] for i in range(1,len(M)+1)]
#chemin surdiagonal -> sous-diagonal
def Sup_to_Sub(d):
    n=len(d)
    if n<=1:
        return d
    a=0
    for k in range(1,n):
        if d[-k-1]==n-k:
            a=n-k
            break
    l=d[a:-1]
    for k in range(n-a-1):
        l[k]-=a+1
    l,m=Sup_to_Sub(l),Sup_to_Sub(d[:a])
    for k in range(n-a-1):
        l[k]+=a+1
    return [1]+m+l
#chemin sous-diagonal -> sur-diagonal
def Sub_to_Sup(d):
    n=len(d)
    if n<=1:
        return d
    a=n
    for k in range(1,n):
        if d[k]==k+1:
            a=k
            break
    l=d[a:]
    for k in range(n-a):
        l[k]-=a
    l,m=Sub_to_Sup(l),Sub_to_Sup(d[1:a])
    for k in range(n-a):
        l[k]+=a
    return m+l+[n]
def Standard_Procedure(G):
    n,G = len(G),[list(g[:]) for g in G]
    for k in range(1,n):
        for i in range(k):
            if G[-i-1][-k-1+i]==G[-i-2][-k+i]:
                for j in range(i+1):
                    G[-i-1+j][-k+i] -= 1
    return G

def Reverse_Standard_Procedure(G):
    n,G = len(G),[g[:] for g in G]
    for k in range(n-1):
        for i in range(1,n-k):
            if G[k+i-1][k]==G[k+i][k]:
                for j in range(n-k-i):
                    G[k+i+j][k+j+1] += 1
    return G

class GT_Tri():
    def __init__(self, T):
        self.main = T
    def __str__(self):
        s,T = '',self.main
        for k in range(len(T)):
            l = []
            for t in T[-k-1]:
                if t==None:
                    l.append('.')
                else:
                    l.append(str(t))
            s += k*" "+" ".join(l)+"\n"
        return s[:-1]
    def __repr__(self):
        return str(self.main)
    def __eq__(self, other):
        if not isinstance(other,GT_Tri):
            return False
        return self.main == other.main
    def __le__(self,other):
        l1,l2 = self.main,other.main
        if len(l1)!=len(l2):
            return False
            raise Exception('les deux triangles doivent avoir les mêmes dimensions')
        for a,b in zip(l1,l2):
            if len(a)!=len(b):
                return False
                raise Exception('les deux triangles doivent avoir les mêmes dimensions')
            for i,j in zip(a,b):
                if i>j:
                    return False
        return True
    def __ge__(self,other):
        l1,l2 = self.main,other.main
        if len(l1)!=len(l2):
            raise Exception('les deux triangles doivent avoir les mêmes dimensions')
        for a,b in zip(l1,l2):
            if len(a)!=len(b):
                raise Exception('les deux triangles doivent avoir les mêmes dimensions')
            for i,j in zip(a,b):
                if i<j:
                    return False
        return True
    def __lt__(self,other):
        return (self<=other) and not (self==other)
    def __gt__(self,other):
        return other < self
    def __hash__(self):
        return repr(self.main).__hash__()
    def draw(self, a = None, size = 70,style = 'num'):
        T = self.main
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        n = len(T)
        if style == None:
            style = 'num'
        if style == 'num':
            for i in range(n):
                for j in range(len(T[i])):
                    if T[i][j]==None:
                        t = '.'
                    else:
                        t = str(T[i][j])
                    ax.text((j-i/2)/n+0.48,i/n*0.9+0.08,t,fontsize=size/n,color='black',ha='center',**tfont)
        elif style == 'til':
            linewidth = size/140
            kmin, kmax = T[-1][0], T[-1][-1]
            for i in range(n-1):
                for j in range(i+1):
                    Draw_Part_Aux(-j,T[i][j]-kmin,i-j,1,ax,linewidth)
                    for k in range(T[i+1][j],T[i][j]):
                        Draw_Part_Aux(-j,k-kmin+1,i-j,0,ax,linewidth)
                    for k in range(T[i][j],T[i+1][j+1]):
                        Draw_Part_Aux(-j-1,k-kmin+1,i-j,2,ax,linewidth)
                for k in range(T[i][-1],kmax):
                    Draw_Part_Aux(-i,k-kmin+1,-1,0,ax,linewidth)
                for k in range(kmin,T[i][0]):
                    Draw_Part_Aux(0,k-kmin+1,i,2,ax,linewidth)
        if a==None:
            plt.show()
    def to_Tab(self):
        T = self.main
        Y,n,z = [],len(T),0
        while z<n and T[-1][z]==0:
            z += 1
        for k in range(n-z):
            y = T[k][0]*[k+1]
            for i in range(1,n-k):
                y += (T[k+i][-k-1]-T[k+i-1][-k-1])*[k+i+1]
            Y.append(y)
        return Young_Tab(Y)
    def Schutz(self):
        T = [t[:] for t in self.main]
        n = len(T)
        for k in range(1,n):
            T = omega(T,n-k)
        return GT_Tri(T)
    def to_ASM(self,width = None):
        G = self.main
        n,m = len(G),G[-1][-1]
        if width:
            m = max(m,width)
        M = [m*[0] for _ in range(n)]
        for i in range(n):
            for j in range(i+1):
                M[i][G[i][j]-1] += 1
        for i in range(n-1,0,-1):
            for j in range(m):
                M[i][j]-=M[i-1][j]
        return ASM(M)
    def to_Bool(self):
        M = self.main
        n = len(M)
        B = [[] for _ in range(n-1)]
        for k in range(1,n):
            a = k
            while a<n and M[a][a]<k+1:
                B[a-1].append(1)
                a += 1
            if a!=n:
                B[a-1].append(0)
            j = 1
            while a<n-1:
                if M[a+1][-j-1]>k:
                    j += 1
                    B[a].append(0)
                else:
                    B[a].append(1)
                a += 1
        return Bool(B)
    def Min(self,other):
        T1,T2 = self.main,other.main
        return GT_Tri([[min(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def Max(self,other):
        T1,T2 = self.main,other.main
        return GT_Tri([[max(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    
class Young_Tab():
    def __init__(self, Y):
        self.main = Y
    def __str__(self):
        s,T = '',self.main
        for k in range(len(T)):
            s += " ".join([str(t) for t in T[-k-1]])+"\n"
        return s[:-1]
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=Young_Tab:
            return False
        return self.main == other.main
    def draw(self, a = None, size = 70,style = None):
        Y = self.main
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        if Y==[]:
            n = 1
        else:
            n = max(len(Y),len(Y[0]))
        for i in range(len(Y)):
            for j in range(len(Y[i])):
                ax.text(j/n,i/n,str(Y[i][j]),fontsize=size/n,color='black',ha='center',**tfont)
        if a==None:
            plt.show()
    def to_Tri(self):
        Y,T = self.main,[]
        for k in range(max(y[-1] for y in Y)):
            t = []
            for i in range(k+1):
                if k-i<len(Y):
                    for j in range(len(Y[k-i])):
                        if Y[k-i][j]>k+1:
                            j -= 1
                            break
                    t.append(j+1)
                else:
                    t.append(0)
            T.append(t)
        return GT_Tri(T)
    def Schutz(self):
        T = self.main
        T2 = [list(t) for t in T]
        n = sum(len(t) for t in T)
        m = max([t[-1] for t in T])+1
        for k in range(n):
            x,y=0,0
            v = T2[0][0]
            while True:
                if y<len(T[x])-1 and T2[x][y+1]>0:
                    if x<len(T)-1 and y<len(T[x+1]) and T2[x+1][y]>0 and T2[x][y+1]>=T2[x+1][y]:
                        T2[x][y] = T2[x+1][y]
                        x+=1
                    else:
                        T2[x][y] = T2[x][y+1]
                        y+=1
                else:
                    if x<len(T)-1 and y<len(T[x+1]) and T2[x+1][y]>0:
                        T2[x][y] = T2[x+1][y]
                        x+=1
                    else:
                        T2[x][y] = -v
                        break
        for i in range(len(T)):
            for j in range(len(T[i])):
                T2[i][j]+=m
        return Young_Tab(T2)
    def Tuple(self):
        return tuple(tuple(l) for l in self.__main)
class ASM():
    def __init__(self, M):
        self.main = M
        self.corner_sum = None
    def __str__(self):
        D = {0:" ",-1: "-",1:"+"}
        s,M = '',self.main
        for l in M:
            s += " ".join([D[t] for t in l])+"\n"
        return s[:-1]
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=ASM:
            return False
        return self.main == other.main
    def __le__(self,other):
        M1,M2 = self.Corner_Sum(), other.Corner_Sum()
        if len(M1)!=len(M2):
            return False
            raise Exception('les deux ASMs doivent avoir les mêmes dimensions')
        for a,b in zip(M1,M2):
            for i,j in zip(a,b):
                if i<j:
                    return False
        return True
    def __ge__(self,other):
        M1,M2 = self.Corner_Sum(), other.Corner_Sum()
        if len(M1)!=len(M2):
            return False
            raise Exception('les deux ASMs doivent avoir les mêmes dimensions')
        for a,b in zip(M1,M2):
            for i,j in zip(a,b):
                if i>j:
                    return False
        return True
    def __lt__(self,other):
        return (self<=other) and not (self==other)
    def __gt__(self,other):
        return other < self
    def Corner_Sum(self):
        if self.corner_sum == None:
            self.corner_sum = Corner_Sum(self.main)
        return self.corner_sum
    def draw(self, a = None, size = 70,style = 'sqr',save = False):
        if style==None:
            style = 'sqr'
        M = self.main
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        n,m = len(M),len(M[0])
        nm = max(n,m)
        if style=='num':
            for i in range(n):
                for j in range(m):
                    if M[i][j]==0:
                        ax.text(j/nm,1-i/nm,str(M[i][j]),fontsize=size/nm,color='black',ha='center',**tfont)
                    elif M[i][j]>0:
                        ax.text(j/nm,1-i/nm,str(M[i][j]),fontsize=size/nm,color='red',fontweight='bold',ha='center',**tfont)
                    else:
                        ax.text(j/nm,1-i/nm,str(M[i][j]),fontsize=size/nm,color='blue',fontweight='bold',ha='center',**tfont)
        if style=='numblack':
            for i in range(n):
                for j in range(m):
                    if M[i][j]==0:
                        ax.text(j/nm,1-i/nm,str(M[i][j]),fontsize=size/nm,color='black',ha='center',**tfont)
                    elif M[i][j]>0:
                        ax.text(j/nm,1-i/nm,str(M[i][j]),fontsize=size/nm,color='black',fontweight='bold',ha='center',**tfont)
                    else:
                        ax.text(j/nm,1-i/nm,str(M[i][j]),fontsize=size/nm,color='black',fontweight='bold',ha='center',**tfont)
        elif style=='sqr':
            for i in range(n):
                for j in range(m):
                    if M[i][j] != None:
                        if isinstance(M[i][j],ColoredInt):
                            color = M[i][j].color
                        else:
                            if M[i][j]==0:
                                color = 'white'
                            elif M[i][j]>0:
                                color = 'tomato'
                            else:
                                color = 'blue'
                        ax.fill([j,j+1,j+1,j],[n-i,n-i,n-i-1,n-i-1],facecolor=color,edgecolor='black',linewidth=size/140)
        if a==None:
            if save:
                title,letters = '',string.ascii_letters
                for _ in range(10):
                    title += random.choice(letters)
                plt.savefig(title+'.svg')
                print(f'File saved as {title}_.svg')
            plt.show()
    
    def to_Tri(self):
        M = self.main
        n,m = len(M),len(M[0])
        l,G = m*[int(0)],[]
        for k in range(n):
            g = []
            for i in range(m):
                l[i] += M[k][i]
                g += l[i]*[i+1]
            G.append(g)
        return GT_Tri(G)
    def Schutz(self):
        M = [m[:] for m in self.main]
        n = len(M)
        for k in range(1,n):
            M = omega_ASM(M,n-k)
        return ASM(M)
    def Tuple(self):
        return tuple(tuple(l) for l in self.main)
    
def Draw_Part_Aux(i,j,k,z,ax,linewidth):
    s3=sqrt(3)/2
    X,Y=[(j-i)*s3],[k-(i+j)/2]
    P=[((j-i)*s3,k-(i+j)/2)]
    if z==0:
        X+=[(j-i+1)*s3,(j-i)*s3,(j-i-1)*s3]
        Y+=[k-(i+j-1)/2,k-(i+j)/2+1,k-(i+j-1)/2]
        ax.fill(X,Y,facecolor='#ffffb3',edgecolor='black',linewidth=linewidth)
    elif z==1:
        X+=[(j-i)*s3,(j-i+1)*s3,(j-i+1)*s3]
        Y+=[k-(i+j)/2-1,k-(i+j+1)/2,k-(i+j-1)/2]
        ax.fill(X,Y,facecolor='#1818ff',edgecolor='black',linewidth=linewidth)
    else:
        X+=[(j-i)*s3,(j-i-1)*s3,(j-i-1)*s3]
        Y+=[k-(i+j)/2-1,k-(i+j+1)/2,k-(i+j-1)/2]
        ax.fill(X,Y,facecolor='#ff6b6b',edgecolor='black',linewidth=linewidth)
        
class Plane_Part():
    """
    sup_part = plane partition containing P, used for drawing
    """
    def __init__(self, P, sup_part = None):
        self.main = P
        self.sup_part = sup_part
    def __str__(self):
        return str(self.main)
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=Plane_Part:
            return False
        return self.main == other.main
    def draw(self, a = None, size = 70,style = None):
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            plt.axis('off')
        else:
            ax=a
        P = self.main
        Q = self.sup_part
        linewidth = size/140
        if Q:
            for i in range(len(Q)):
                for j in range(len(Q[i])):
                    if i<len(P) and j<len(P[i]):
                        Draw_Part_Aux(i,j,P[i][j]-1,0,ax,linewidth)
                        if j==len(P[i])-1:
                            t=0
                        else:
                            t=P[i][j+1]
                        for k in range(t,P[i][j]):
                            Draw_Part_Aux(i,j,k,1,ax,linewidth)
                        if i==len(P)-1 or j>=len(P[i+1]):
                            t=0
                        else:
                            t=P[i+1][j]
                        for k in range(t,P[i][j]):
                            Draw_Part_Aux(i,j,k,2,ax,linewidth)
                    else:
                        Draw_Part_Aux(i,j,-1,0,ax,linewidth)
                if i<len(P):
                    t = P[i][0]
                else:
                    t = 0
                for k in range(t,Q[i][0]):
                    Draw_Part_Aux(i,-1,k,1,ax,linewidth)
            for j in range(len(Q[0])):
                if len(P)>0 and j<len(P[0]):
                    t = P[0][j]
                else:
                    t = 0
                for k in range(t,Q[0][j]):
                    Draw_Part_Aux(-1,j,k,2,ax,linewidth)
        for i in range(len(P)):
            for j in range(len(P[i])):
                Draw_Part_Aux(i,j,P[i][j]-1,0,ax,linewidth)
                if j==len(P[i])-1:
                    t=0
                else:
                    t=P[i][j+1]
                for k in range(t,P[i][j]):
                    Draw_Part_Aux(i,j,k,1,ax,linewidth)
                if i==len(P)-1 or j>=len(P[i+1]):
                    t=0
                else:
                    t=P[i+1][j]
                for k in range(t,P[i][j]):
                    Draw_Part_Aux(i,j,k,2,ax,linewidth)
        if a==None:
            plt.show()

class DPP():
    def __init__(self, P,order):
        self.main = P
        self.order = order
    def __str__(self):
        s,T = '',self.main
        for k in range(len(T)):
            s += 2*k*" "+" ".join([str(t) for t in T[k]])+"\n"
        return s[:-1]
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=DPP:
            return False
        return self.main == other.main
    def __le__(self,other):
        if self.order!=other.order:
            raise Exception('les DPPs doivent être du même ordre')
        l1,l2 = self.main,other.main
        if len(l1)>len(l2):
            return False
        for a,b in zip(l1,l2):
            if len(a)>len(b):
                return False
            for i,j in zip(a,b):
                if i>j:
                    return False
        return True
    def __ge__(self,other):
        return other <= self
    def __lt__(self,other):
        return (self<=other) and not (self==other)
    def __gt__(self,other):
        return other < self
    def draw(self, a = None, size = 70,style = 'num'):
        P,n = self.main,self.order-1
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        if style == None:
            style = 'num'
        if style == 'num':
            for i in range(len(P)):
                for j in range(len(P[i])):
                    if P[i][j]==None:
                        t = ''
                    else:
                        t = str(P[i][j])
                    ax.text((j+i)/n,-1.1*i/n+1,t,fontsize=size/n,color='black',ha='center',**tfont)
        elif style == 'til':
            linewidth = size/140
            P = [l+ (l[0]-len(l)-1)*[0] for l in P]
            for i in range(len(P)):
                for j,k in enumerate(P[i]):
                    Draw_Part_Aux(0,j,k,0,ax,linewidth)
                    Draw_Part_Aux(k,0,j+2,2,ax,linewidth)
                    Draw_Part_Aux(j+2,k,2,1,ax,linewidth)
                for j in range(1,len(P[i])):
                    for k in range(P[i][j],P[i][j-1]):
                        Draw_Part_Aux(0,j-1,k+1,1,ax,linewidth)
                        Draw_Part_Aux(k+1,0,j+1,0,ax,linewidth)
                        Draw_Part_Aux(j+1,k+1,2,2,ax,linewidth)
                    k0 = P[i+1][j-1]+1 if i+1<len(P) and j<=len(P[i+1]) else 0
                    for k in range(k0,P[i][j]):
                        Draw_Part_Aux(0,j,k+1,2,ax,linewidth)
                        Draw_Part_Aux(k+1,0,j+2,1,ax,linewidth)
                        Draw_Part_Aux(j+2,k+1,2,0,ax,linewidth)
                for k in range(P[i][-1]):
                    Draw_Part_Aux(0,len(P[i])-1,k+1,1,ax,linewidth)
                    Draw_Part_Aux(k+1,0,len(P[i])+1,0,ax,linewidth)
                    Draw_Part_Aux(len(P[i])+1,k+1,2,2,ax,linewidth)
            l0 = list(P[0]) if len(P)>0 else []
            l0 += (n-len(l0))*[-1]
            for j in range(n):
                for k in range(l0[j],n+1):
                    Draw_Part_Aux(-1,j,k+1,2,ax,linewidth)
                    Draw_Part_Aux(k+1,-1,j+2,1,ax,linewidth)
                    Draw_Part_Aux(j+2,k+1,1,0,ax,linewidth)
        if a==None:
            plt.show()
    def dual(self):
        N,P = self.order,self.main
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
    def to_Tri(self):
        P = self.main
        n = self.order
        T = []
        for k in range(n):
            t = []
            for i in range(k+1):
                if k-i<len(P):
                    for j in range(len(P[k-i])):
                        if P[k-i][j]<n-k:
                            j -= 1
                            break
                    t.append(j+2)
                else:
                    t.append(1)
            T.append(t)
        return GT_Tri(T)

class Bool():
    def __init__(self, T):
        self.main = T
    def __str__(self):
        s,T = '',self.main
        for k in range(len(T)):
            s += k*" "+" ".join([str(t) for t in T[-k-1]])+"\n"
        return s[:-1]
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=Bool:
            return False
        return self.main == other.main
    def __le__(self,other):
        l1,l2 = self.main,other.main
        if len(l1)!=len(l2):
            return False
            raise Exception('les deux triangles doivent avoir les mêmes dimensions')
        for a,b in zip(l1,l2):
            if len(a)!=len(b):
                return False
                raise Exception('les deux triangles doivent avoir les mêmes dimensions')
            for i,j in zip(a,b):
                if i>j:
                    return False
        return True
    def __ge__(self,other):
        l1,l2 = self.main,other.main
        if len(l1)!=len(l2):
            raise Exception('les deux triangles doivent avoir les mêmes dimensions')
        for a,b in zip(l1,l2):
            if len(a)!=len(b):
                raise Exception('les deux triangles doivent avoir les mêmes dimensions')
            for i,j in zip(a,b):
                if i<j:
                    return False
        return True
    def draw(self, a = None, size = 70,style = None):
        T = self.main
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        n = len(T)
        for i in range(n):
            for j in range(len(T[i])):
                if T[i][j]==None:
                    t = ''
                else:
                    t = str(T[i][j])
                ax.text((j-i/2)/n+0.48,i/n*0.9+0.08,t,fontsize=size/n,color='black',ha='center',**tfont)
        if a==None:
            plt.show()
    def to_Magog(self):
        B = [l[:] for l in self.__main]
        n = len(B)
        for k in range(n):
            for i in range(k):
                B[k][i] += B[k-1][i]
        for k in range(n):
            for i in range(k):
                B[k][-i-1] += 1-B[k][-i-2]
            B[k][0] += 1
        M = [[1]]
        for k in range(n):
            l = []
            for i in range(k+1):
                l += B[k][i]*[i+1]
            if len(l)==k+1:
                l.append(k+2)
            M.append(l)
        return GT_Tri(M)

class List():
    def __init__(self, l):
        self.main = l
    def __str__(self):
        return " ".join([str(k) for k in self.main])
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        return self.main == other.main
    def __le__(self,other):
        l1,l2 = self.main,other.main
        if len(l1)!=len(l2):
            raise Exception('les deux listes doivent avoir les mêmes dimensions')
        for i,j in zip(l1,l2):
            if i>j:
                return False
        return True
    def __ge__(self,other):
        return other <= self
    def __lt__(self,other):
        return (self<=other) and not (self==other)
    def __gt__(self,other):
        return other < self
    def draw(self, a = None, size = 70,style = None):
        l = self.main
        n = len(l)
        ASM([[int(i>=n-l[j]) for j in range(n)] for i in range(n)]).draw(a,size,style)
class Subpath(List):
    def Down_Covers(self):
        P = self.main
        n,C = len(P),[]
        for k in range(1,n):
            if P[k]>P[k-1]:
                P2 = P.copy()
                P2[k] -= 1
                C.append(Subpath(P2))
        return C
    def Up_Covers(self):
        P = self.main
        n,C = len(P),[]
        for k in range(1,n-1):
            if P[k]<min(P[k+1],k+1):
                P2 = P.copy()
                P2[k] += 1
                C.append(Subpath(P2))
        if P[-1]<n:
            P2 = P.copy()
            P2[-1] += 1
            C.append(Subpath(P2))
        return C
    def Min(self,other):
        return Subpath([min(i,j) for i,j in zip(self.main,other.main)])
    def Max(self,other):
        return Subpath([max(i,j) for i,j in zip(self.main,other.main)])
class Parking(List):
    def Down_Covers(self):
        P = self.main
        n,C = len(P),[]
        for k in range(1,n):
            if P[k]>1:
                P2 = P.copy()
                P2[k] -= 1
                C.append(Parking(P2))
        return C
    def Up_Covers(self):
        P = self.main
        n,C = len(P),[]
        for k in range(1,n):
            if P[k]<k+1:
                P2 = P.copy()
                P2[k] += 1
                C.append(Parking(P2))
        return C
    def Min(self,other):
        return Parking([min(i,j) for i,j in zip(self.main,other.main)])
    def Max(self,other):
        return Parking([max(i,j) for i,j in zip(self.main,other.main)])
class Path(List):
    def __init__(self, l, m = None):
        self.main = l
        self.max = m
    def Down_Covers(self):
        P = self.main
        n,C = len(P),[]
        if P[0]>0:
            P2 = P.copy()
            P2[0] -= 1
            C.append(Path(P2))
        for k in range(1,n):
            if P[k]>P[k-1]:
                P2 = P.copy()
                P2[k] -= 1
                C.append(Path(P2))
        return C
    def Up_Covers(self):
        P = self.main
        n,C = len(P),[]
        for k in range(n-1):
            if P[k]<P[k+1]:
                P2 = P.copy()
                P2[k] += 1
                C.append(Path(P2,self.max))
        if P[-1]<self.max:
            P2 = P.copy()
            P2[-1] += 1
            C.append(Path(P2,self.max))
        return C
    def Min(self,other):
        return Path([min(i,j) for i,j in zip(self.main,other.main)])
    def Max(self,other):
        return Path([max(i,j) for i,j in zip(self.main,other.main)])

#P = couple de listes de 0 et de 1
class NonIntersectingPath():
    def __init__(self, P):
        self.main = P
    def main(self):
        return self.main
    def __str__(self):
        return str(self.main)
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=NonIntersectingPath:
            return False
        return self.main == other.main
    def __le__(self,other):
        (l1,l2),(m1,m2) = self.main,other.main
        s1,s2,t1,t2 = 0,0,0,0
        for i,j in zip(l1,m1):
            s1 += i
            t1 += j
            if s1>t1:
                return False
        for i,j in zip(l2,m2):
            s2 += i
            t2 += j
            if s2>t2:
                return False
        return True
    def __ge__(self,other):
        return other <= self
    def __lt__(self,other):
        return not (other <= self)
    def __gt__(self,other):
        return (self <= other)
def NonIntersectingPaths(n):
    L = [[([],[])]]
    for k in range(n):
        L2 = []
        for l1,l2 in L[k]:
            L2.append(([0]+l1,[0]+l2))
            L2.append(([1]+l1,[1]+l2))
        for i in range(k):
            for l1,l2 in L[i]:
                for m1,m2 in L[k-i-1]:
                    L2.append(([0]+l1+[1]+m1,[1]+l2+[0]+m2))
        L.append(L2)
    return [NonIntersectingPath(l) for l in L[-1]]

#codé sous forme de liste croissante de trois entiers distincts entre 0 et n
class Bigrassmannian():
    def __init__(self, P,typ = 'gog',size = None):
        self.main = P
        self.type = typ
        self.size = size
    def __str__(self):
        return str(self.main)
    def __repr__(self):
        return str(self.main)
    def __hash__(self):
        return str(self.main).__hash__()
    def __eq__(self, other):
        if type(other)!=Bigrassmannian:
            return False
        return self.main == other.main and self.size == other.size
    def __le__(self,other):
        P,Q = self.main,other.main
        if self.type == 'gog':
            return P[0]>=Q[0] and P[2]<=Q[2] and Q[0]-P[0]<=Q[1]-P[1]<=Q[2]-P[2]
        elif self.type == 'magog':
#            return P[1]-P[0]-1<=Q[1]-Q[0]-1 and P[2]<=Q[2] and 0<=Q[1]-P[1]<=Q[2]-P[2]
            return P[0]<=Q[0] and 0<=Q[1]-P[1]<=Q[2]-P[2]
        elif self.type == 'gapless':
            return 0<=Q[2]-P[2] and 0<=Q[0]-P[0]<=Q[1]-P[1]
        elif self.type == 'catalan':
            return max(0,Q[0]-P[0])<=min(Q[1]-P[1],Q[2]-P[2])
        elif self.type == 'GTtri':
            return 0<=Q[0]-P[0]<=Q[2]-P[2]<=Q[1]-P[1]
    def __ge__(self,other):
        return other <= self
    def __lt__(self,other):
        return (self<=other) and not (self==other)
    def __gt__(self,other):
        return other < self
    def center(self):
        i,j,k = self.main
        if self.type == 'gog':
            return self.size - j,k-j+i
        if self.type == 'magog':
            return j-i-1,k-j
    def to_Tri(self,mode = 'join'):
        if self.type == 'gog':
            i,j,k = self.main
            G,n = [],self.size
            if mode=='meet':
                i,k = j-i-1,j+n-k+1
                for a in range(i):
                    G.append(list(range(n-a,n+1)))
                if i>0:
                    for a in range(j-i):
                        G.append(list(range(n-i-k+j-a,n-i-k+j+1))+G[i-1].copy())
                    for a in range(k-j):
                        G.append(G[j-1][:j-i]+list(range(n-a-i,n-i+1))+G[i-1].copy())
                else:
                    for a in range(j):
                        G.append(list(range(n-a-k+j,n-k+j+1)))
                    for a in range(k-j):
                        G.append(G[j-1][:j]+list(range(n-a,n+1)))
                for a in range(k,n):
                    G.append(list(range(n-a,n+1)))
            elif mode=='join':
                for a in range(i):
                    G.append(list(range(1,a+2)))
                if i>0:
                    for a in range(j-i):
                        G.append(G[i-1].copy()+list(range(i+k-j+1,i+k-j+a+2)))
                    for a in range(k-j):
                        G.append(G[i-1].copy()+list(range(i+1,a+i+2))+G[j-1][i:])
                else:
                    for a in range(j):
                        G.append(list(range(k-j+1,k-j+a+2)))
                    for a in range(k-j):
                        G.append(list(range(1,a+2))+G[j-1])
                for a in range(k,n):
                    G.append(list(range(1,a+2)))
            else:
                raise Exception("mode must be 'meet' or 'join'")
            return Gog(G)
        elif self.type == 'magog':
            i,j,k = self.main
            M,n = [],self.size
            if mode=='join':
                for a in range(n-j):
                    M.append((a+1)*[1])
                for a in range(i+1):
                    M.append((n-j)*[1]+(a+1)*[k-j+1])
                for a in range(j-i-1):
                    M.append((n-j+a+1)*[1]+(i+1)*[k-j+1])
            elif mode=='meet':
                for a in range(i-j+k):
                    M.append(list(range(1,a+2)))
                for a in range(n-k):
                    M.append(list(range(1,k-j+1))+(a+1)*[k-j]+list(range(k-j+a+2,i-j+k+a+2)))
                for a in range(j-i):
                    M.append(list(range(1,k-j+1))+(n-k+1)*[k-j]+list(range(n-j+2,n+i-j+a+2)))
            return Magog(M)
        elif self.type == 'gapless':
            i,j,k = self.main
            T,n = [],self.size
            if mode=='join':
                for a in range(n):
                    t = []
                    for b in range(a+1):
                        if a-b<=i:
                            if a<i+k-j:
                                c = j-i+1+b
                            else:
                                c = b-a+k
                            t.append(max(c,b+1))
                        else:
                            t.append(b+1)
                    T.append(t)
            if mode=='meet':
                for a in range(n):
                    t = []
                    for b in range(a+1):
                        if a-b>=i:
                            if a<i+k-j:
                                c = k-a+b-1
                            else:
                                c = b+j-i
                            t.append(min(c,n-a+b))
                        else:
                            t.append(n-a+b)
                    T.append(t)
            return Gog(T)
        elif self.type == 'catalan':
            i,j,k = self.main
            T,n = [],self.size
            if mode=='join':
                for a in range(n):
                    t = []
                    nl = k-1-j+i
                    for b in range(a+1):
                        if b>i and b-a>1+j-k:
                            c = j+1
                        elif b<i and b-a<1+j-k:
                            c = k-i+2*b-a
                        elif a>nl:
                            c = k-a+b
                        else:
                            c = j+1-i+b
                        t.append(max(c,b+1))
                    T.append(t)
            elif mode=='meet':
                for a in range(n):
                    t = []
                    nl = k-1-j+i
                    for b in range(a+1):
                        if b>i and b-a>1+j-k:
                            c = k-i+2*b-1-a
                        elif b<i and b-a<1+j-k:
                            c = j
                        elif a>nl:
                            c = j-i+b
                        else:
                            c = k-1-a+b
                        t.append(min(c,n-a+b))
                    T.append(t)                    
            return CatalanTriangle(T)
        elif self.type == 'GTtri':
            i,j,k = self.main
            T,n = [],self.size
            if mode=='join':
                for a in range(n):
                    t = []
                    nl = k-1-j+i
                    for b in range(a+1):
                        if b>=k-j-1 and a-b<=i:
                            c = k-i
                        else:
                            c = 0
                        t.append(max(c,b+1))
                    T.append(t)
            elif mode=='meet':
                for a in range(n):
                    t = []
                    nl = k-1-j+i
                    for b in range(a+1):
                        if b<=k-j-1 and a-b>=i:
                            c = k-i-1
                        else:
                            c = n
                        t.append(min(c,n-a+b))
                    T.append(t)                    
            return GT_Tri(T)
        
def Bigrassmannians(n,typ = 'gog'):
    return [Bigrassmannian((i,j,k),typ,n) for i in range(n-1) for j in range(i+1,n) for k in range(j+1,n+1)]

class Gog(GT_Tri):
    def __init__(self, T):
        if isinstance(T,GT_Tri):
            self.main = T.main
        else:
            self.main = T
    def upper_covers(self):
        T = self.main
        n,C = len(T),[]
        for k in range(n-1):
            for i in range(k):
                if min(T[k-1][i],min(T[k+1][i+1],T[k][i+1]-1))>T[k][i]:
                    T2 = [t.copy() for t in T]
                    T2[k][i] += 1
                    C.append(Gog(T2))
            if T[k+1][-1]>T[k][-1]:
                T2 = [t.copy() for t in T]
                T2[k][-1] += 1
                C.append(Gog(T2))
        return C
    def lower_covers(self):
        T = self.main
        n,C = len(T),[]
        for k in range(n-1):
            if T[k][0]>T[k+1][0]:
                T2 = [t.copy() for t in T]
                T2[k][0] -= 1
                C.append(Gog(T2))
            for i in range(1,k+1):
                if T[k][i]>max(T[k-1][i-1],max(T[k+1][i],T[k][i-1]+1)):
                    T2 = [t.copy() for t in T]
                    T2[k][i] -= 1
                    C.append(Gog(T2))
        return C
    def Min(self,other):
        T1,T2 = self.main,other.main
        return Gog([[min(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def Max(self,other):
        T1,T2 = self.main,other.main
        return Gog([[max(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def canonical_joinands(self, form = 'elem'):
        n = len(self.main)
        if form not in ['elem', 'irr']:
            raise Exception("form should be 'elem' or 'irr'")
        l = []
        for b in Bigrassmannians(n):
            B = b.to_Tri('join')
            if self >= B:
                if form == 'elem':
                    l.append(B)
                elif form == 'irr':
                    l.append(b)
        return [x for x in l if all(not x < y for y in l)]
    def canonical_meetands(self, form = 'elem'):
        n = len(self.main)
        if form not in ['elem', 'irr']:
            raise Exception("form should be 'elem' or 'irr'")
        l = []
        for b in Bigrassmannians(n):
            B = b.to_Tri('meet')
            if self <= B:
                if form == 'elem':
                    l.append(B)
                elif form == 'irr':
                    l.append(b)
        return [x for x in l if all(not x > y for y in l)]
    
class Magog(GT_Tri):
    def __init__(self, T):
        if isinstance(T,GT_Tri):
            self.main = T.main
        else:
            self.main = T
    def upper_covers(self):
        T = self.main
        n,C = len(T),[]
        for k in range(1,n-1):
            for i in range(1,k):
                if min(T[k+1][i+1],min(T[k-1][i],i+1))>T[k][i]:
                    T2 = [t.copy() for t in T]
                    T2[k][i] += 1
                    C.append(Magog(T2))
            if min(k+1,T[k+1][-1])>T[k][-1]:
                T2 = [t.copy() for t in T]
                T2[k][-1] += 1
                C.append(Magog(T2))
        for k in range(n-1):
            if min(T[-2][k],k+1)>T[-1][k]:
                T2 = [t.copy() for t in T]
                T2[-1][k] += 1
                C.append(Magog(T2))
        if n>T[-1][-1]:
            T2 = [t.copy() for t in T]
            T2[-1][-1] += 1
            C.append(Magog(T2))
        return C
    def lower_covers(self):
        T = self.main
        n,C = len(T),[]
        for k in range(n-1):
            for i in range(k):
                if T[k][i+1]>max(T[k-1][i],T[k+1][i+1]):
                    T2 = [t.copy() for t in T]
                    T2[k][i+1] -= 1
                    C.append(Magog(T2))
        for k in range(1,n):
            if T[-1][k]>T[-2][k-1]:
                T2 = [t.copy() for t in T]
                T2[-1][k] -= 1
                C.append(Magog(T2))
        return C
    def Min(self,other):
        T1,T2 = self.main,other.main
        return Magog([[min(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def Max(self,other):
        T1,T2 = self.main,other.main
        return Magog([[max(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def to_TSSCPP(self):
        M = self.main
        P,n = [],len(M)
        for i in range(1,n):
            p = []
            for j in range(1,n):
                a = M[n-max(i,j)+min(i,j)-1][n-max(i,j)]
                if a>1:
                    p.append(a-1)
                else:
                    break
            if p==[]:
                break
            else:
                P.append(p)
        T = [n*[2*n]+n*[n] for _ in range(n)]+[n*[n] for _ in range(n)]
        for i in range(len(P)):
            T[n+i] += P[i]
            for j in range(len(P[i])):
                T[n-i-1][n-j-1] -= P[i][j]
        P = Rotate(P)
        for i in range(len(P)):
            for j in range(len(P[i])):
                T[n+i][j] += P[i][j]
                T[n-i-1][2*n-j-1] -= P[i][j]
        P = Rotate(P)
        for i in range(len(P)):
            for j in range(len(P[i])):
                T[i][n+j] += P[i][j]
                T[2*n-i-1][n-j-1] -= P[i][j]
        return Plane_Part(T,[(2*n)*[2*n] for _ in range(2*n)])

class AST(GT_Tri):
    def __init__(self, T):
        if isinstance(T,GT_Tri):
            self.main = T.main
        else:
            self.main = T
    def Min(self,other):
        T1,T2 = self.main,other.main
        return AST([[min(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def Max(self,other):
        T1,T2 = self.main,other.main
        return AST([[max(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def to_ASM(self):
        G = self.main
        n = len(G)
        m = 2*n-1
        M = [m*[int(0)] for _ in range(n)]
        for i in range(n):
            for j in range(n-i):
                M[i][G[-i-1][j]-1] += int(1)
        for i in range(1,n):
            for j in range(m):
                M[i-1][j] -= M[i][j]
        for i in range(n-1):
            for j in range(n-i-1):
                M[i][j] = None
                M[i][-j-1] = None
        return ASM(M)
    def upper_covers(self):
        cov = []
        T = self.main
        n = len(T)
        for row in range(n):
            for col in range(row + 1):
                if row < n-1:
                    upbound = T[row + 1][col + 1]
                    if row-col != 0:
                        upbound = min(upbound, T[row][col + 1] - 1)
                        upbound = min(upbound, T[row - 1][col])
                    if T[row][col] < upbound:
                        nv = T[row][col] + 1
                        if (row-col != 0 and nv == T[row - 1][col]) or row+1 <= nv <= 2*n - row - 1:
                            if T[row][col] != T[row + 1][col] or row+2 <= T[row + 1][col] <= 2*n - row-2:
                                T2 = [l.copy() for l in T]
                                T2[row][col] += 1
                                cov.append(AST(T2))
                else:
                    v = T[row][col]
                    if v < n-1:
                        if v in T[v-1] and (v == 1 or v not in T[v-2]) and all(v+1 not in l for l in T[v-1:v+1]):
                            T2 = [l.copy() for l in T]
                            for i in range(v-1, n):
                                k = T[i].index(v)
                                T2[i][k] = v+1
                            cov.append(AST(T2))
                    elif v == n-1:
                        if v in T[n-2] and (n==2 or v not in T[n-3]) and n not in T[n-2] and n+1 not in T[n-2]:
                            T2 = [l.copy() for l in T]
                            k = T[n-2].index(v)
                            T2[n-2][k] = v+2
                            T2[row][col] = n
                            T2[row][col+1] = n+1
                            cov.append(AST(T2))
                    elif n+1 <= v < 2*n-1:
                        if v in T[2*n-v-2] and v+1 not in T[2*n-v-2]:
                            T2 = [l.copy() for l in T]
                            for i in range(2*n-v-2, n):
                                k = T[i].index(v)
                                T2[i][k] = v+1
                            cov.append(AST(T2))   
        return cov

class CatalanTriangle(GT_Tri):
    def __init__(self, T):
        if isinstance(T,GT_Tri):
            self.main = T.main
        else:
            self.main = T
    def Min(self,other):
        T1,T2 = self.main,other.main
        return CatalanTriangle([[min(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def Max(self,other):
        T1,T2 = self.main,other.main
        return CatalanTriangle([[max(i,j) for i,j in zip(l1,l2)] for l1,l2 in zip(T1,T2)])
    def to_TriangleSubposet(self):
        r"""
        Returns a list of bigrassmannians whose subposet is isomorphic to a triangle poset.
        """
        lb,n = [],len(self.main)
        for i,l in enumerate(self.main):
            for j,k in enumerate(l):
                lb.append(Bigrassmannian((j,n+1-i+2*j-k ,n+1-i+j),size = n+1))
        return lb
    def Sym(self):
        n = len(self.main)
        return CatalanTriangle([[2*j-i+k for j,k in enumerate(reversed(l))] for i,l in enumerate(self.main)])
    def Dual(self):
        n = len(self.main)
        return CatalanTriangle([[n+1-k for j,k in enumerate(reversed(l))] for i,l in enumerate(self.main)])
    def to_Matrix(self):
        C = self.main
        n = len(C)-1
        M = [n*[0] for _ in range(n)]
        for i in range(n):
            for j in range(i+1):
                if i==n-1 or C[i][j]==C[i+2][j+1]:
                    if C[i][j]==C[i+1][j]!=C[i+1][j+1]:
                        M[j][j-i+n-1] = -1
                    elif C[i][j]==C[i+1][j+1]!=C[i+1][j]:
                        M[j][j-i+n-1] = 1
        return ASM(M)
    def to_PipeDream(self):
        r"""
        Needs Pipe Dreams.sage
        """
        T = []
        C = self.main
        n = len(C)-1
        c = ['#ff0000' if C[-1][-1]==C[-2][-1] else '#00ffff']
        for i in range(1,n):
            t = []
            for j in range(i):
                if C[-i+j][n-i]==C[-i+j-2][n-i-1] and C[-i+j-1][n-i]==C[-i+j-1][n-i-1]+1:
                    t.append(1)
                else:
                    t.append(0)
            T.append(t)
            if C[-i-1][-1]==C[-i-2][-1]:
                c.append('#ff0000')
            else:
                c.append('#00ffff')
        return PipeDream(T,c)
        
class ColoredInt(int):
    def __new__(cls, n, color):
        if isinstance(n,ColoredInt):
            x = int.__new__(cls,n)
            x.color = n.color
        else:
            x = int.__new__(cls,n)
            x.color = color
        return x

class ColoredGog():
    def __init__(self, T):
        self.main = T
    def __str__(self):
        s,T = '',self.main
        for k in range(len(T)):
            s += k*" "+" ".join([str(t) for t in T[-k-1]])+"\n"
        return s[:-1]
    def __repr__(self):
        return str(self.main)
    def __eq__(self, other):
        if type(other) == ColoredGog:
            return self.main == other.main
        return False
    def __le__(self,other):
        l1,l2 = self.main,other.main
        if len(l1)!=len(l2):
            raise Exception('les deux triangles doivent avoir les mêmes dimensions')
        for a,b in zip(l1,l2):
            if len(a)!=len(b):
                raise Exception('les deux triangles doivent avoir les mêmes dimensions')
            for i,j in zip(a,b):
                if i[0]>j[0]:
                    return False
                if i[0]==j[0] and i[1]!=None and j[1]!=None and i[1]>j[1]:
                    return False
        return True
    def __ge__(self,other):
        l1,l2 = self.main,other.main
        if len(l1)!=len(l2):
            raise Exception('les deux triangles doivent avoir les mêmes dimensions')
        for a,b in zip(l1,l2):
            if len(a)!=len(b):
                raise Exception('les deux triangles doivent avoir les mêmes dimensions')
            for i,j in zip(a,b):
                if i[0]<j[0]:
                    return False
                if i[0]==j[0] and i[1]!=None and j[1]!=None and i[1]<j[1]:
                    return False
        return True
    def __lt__(self,other):
        return (self<=other) and not (self==other)
    def __gt__(self,other):
        return other < self
    def __hash__(self):
        return repr(self.main).__hash__()
    def draw(self, a = None, size = 70,style = None):
        T = self.main
        if a==None:
            fig, ax = plt.subplots(1,1, figsize=(5,5),subplot_kw={'aspect': 'equal'})
            ax.axis('off')
        else:
            ax = a
        n = len(T)
        for i in range(n):
            for j in range(len(T[i])):
                t = str(T[i][j][0])
                if T[i][j][1]==0:
                    c = 'blue'
                elif T[i][j][1]==1:
                    c = 'red'
                elif T[i][j][1]==2:
                    c = 'green'
                else:
                    c = 'black'
                ax.text((j-i/2)/n+0.48,i/n*0.9+0.08,t,fontsize=size/n,color=c,ha='center',**tfont)
        if a==None:
            plt.show()
