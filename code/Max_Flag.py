#!/usr/bin/env python
# coding: utf-8

# # Header

# In[19]:


get_ipython().run_line_magic('run', 'Classes_Weak_Order_Quotients.ipynb')
get_ipython().run_line_magic('run', 'LudoSemiDIst.ipynb')
camb_4 = [([],1,3,[2]),([3],2,4,[])]
def bool_lattice(n):
    return [([i],i-1,i+1,[])for i in range(2,n)] + [([],i-1,i+1,[i])for i in range(2,n)]


# In[ ]:


permu_order_poset = Poset({0:[1,2],1:[3],2:[3]})


# # Functions

# In[2]:


def is_maxflag(L, ret = False):
    P = L.join_irreducibles_poset()
    E, E2 = [], set()
    for x in L:
        m = P.subposet([j for j in P if L.is_lequal(j,x)]).maximal_elements()
        if len(m) == 2:
            E.append(m)
        elif len(m) >= 2:
            E2.add(tuple(sorted(m)))
    G = Graph(E)
    def is_simplicial():
        for m in E2:
            for k, i in enumerate(m):
                for j in m[k+1:]:
                    if not G.has_edge(i,j):
                        return False
        return True
    def is_flag():
        for c in sage.graphs.cliquer.all_cliques(G, 3):
            cs = tuple(sorted(c))
            if cs not in E2:
                return False
        return True
    s = ''
    if not is_simplicial():
        s += 'not simplicial  '
    if not is_flag():
        s += 'not flag  '
    if s != '':
        return s
    if ret:
        return E
    return True


# In[3]:


def forcingOrder(n,essential=True):
    all_arcs = SingleArc.allArcs(n,essential=essential)
    vertices = []
    rels = []
    for i in range(len(all_arcs)):
        #print(type(all_arcs[i]))
        vertices.append(all_arcs[i])
        for j in range(len(all_arcs)):
            if i == j:
                pass
            else:
                if all_arcs[i].isSubarc(all_arcs[j]):
                    rels.append([all_arcs[i],all_arcs[j]])
    forcing_order = Poset((vertices,rels))
    return forcing_order



# In[4]:


def allPermuTreesClasses(n):
    s = [i for i in range(n-2)]
    res = [[]]
    words = []
    for ele in s:
        ns = [j+[ele] for j in res]
        res.extend(ns)

    for not_2 in res:
        #print("not_2:",not_2)
        nb_not2 = len(not_2)
        if nb_not2 == 0:
            words.append([4]+[2 for j in range(n-2)] +[4])
        else:
            xs = [[]]
            for ele in not_2:
                ns = [j+[ele] for j in xs]
                xs.extend(ns)
            for positions in xs:
                nb_xs = len(positions)
                if nb_xs == 0:
                    word=[]
                    not_2idx = 0
                    for idx in range(n-2):
                        if not_2idx < nb_not2 and not_2[not_2idx] == idx:
                            not_2idx+=1
                            word.append(1)
                        else:
                            word.append(2)
                    #print("word:",word)
                    words.append([4]+word+[4])
                else:
                    #print("pos:",positions)         
                    word = []
                    not_2idx = 0
                    pos_idx = 0
                    for idx in range(n-2):
                        if not_2idx < nb_not2 and not_2[not_2idx] == idx:
                            not_2idx+=1
                            if pos_idx < nb_xs and positions[pos_idx] == idx:
                                pos_idx+=1
                                word.append(4)
                            else:
                                word.append(1)
                        else:
                            word.append(2)
                    #print("word",word)
                    words.append([4]+word+[4])
    return words




# In[5]:


def allPermutree(n):
    decorations = [0,1,2,3]
    list_of_permutrees = [[]]
    for i in range(1,n-1):
        list_of_permutrees = [word +[j] for word in list_of_permutrees for j in decorations]
    return list_of_permutrees


# In[10]:


def permuQuotient(decoration,n):
    arc_list = []
    for i in range(n-2):
        #print(i)
        #The decorations of the first and last node have no bearing on the resulting quotient
        node = decoration[i]
        if node in [1,3]:
            arc_list.append(([],i+1,i+3,[i+2]))
        if node in [2,3]:
            arc_list.append(([i+2],i+1,i+3,[]))
    #print(arc_list)
    #print(n)
    return WOQuotient.fromArcList(arc_list,n)


# In[11]:


#returns if permu_A is lesser (or equal) than permu_B
def permu_order(permu_A,permu_B):
    return all (permu_order_poset.is_lequal(i,j) for i,j in zip(permu_A,permu_B))


# In[12]:


def ludoWord(decoration):
    dico = {0:"(|)",1:"(V)",2:"(^)",3:"(X)"}
    return "".join([dico[i] for i in decoration])



# In[13]:


def biWordToMaxFlagPermu(bin_word):
    state = 0
    permu = []
    for char in bin_word:
        #print("char",char,"state",state)
        if state == 0:
            if char == 0:
                permu.append(2)
            else:
                permu.append(0)
                state = 1
        else:
            if char == 0:
                permu.append(3)
                state = 0
            else:
                permu.append(1)

    if state == 0:
        permu.append(0)
    else:
        permu.append(1)
    return permu


# In[14]:


def allMaxFlags(n):
    max_flags = [[0],[1],[2],[3]]
    for i in range(n-3):
        max_flags2 = []
        for elt in max_flags:
            if elt[-1] == 0 or elt[-1] == 1:
                max_flags2.append(elt+[1])
                max_flags2.append(elt+[3])
            else:
                max_flags2.append(elt+[0])
                max_flags2.append(elt+[1])
                max_flags2.append(elt+[2])
                max_flags2.append(elt+[3])
        max_flags = max_flags2

    return max_flags


"""
# # Tests

# In[15]:


n=7


# In[16]:


for t in Tuples([0,1],n-3):
    deco = biWordToMaxFlagPermu(t)
    quo_lattice = LatticePoset(permuQuotient(deco,n).poset())
    JI_poset = quo_lattice.join_irreducibles_poset()
    print(t)
    #JI_poset.show()
    print(JI_poset.antichains().cardinality())
    print(len(JI_poset.list()))
    print()



# In[17]:


list_of_length_max_flag = []
for t in Tuples([0,1],n-3):
    deco = biWordToMaxFlagPermu(t)
    list_of_length_max_flag.append(len(permuQuotient(deco,n).minElements()))
    print(t,end=", ")
    print(list_of_length_max_flag[-1],end=", ")
    print(Integer(list_of_length_max_flag[-1]).factor(),end="")
    print()


# In[18]:


i=0
for deco in allMaxFlags(n):
    card = len(permuQuotient(deco,n).minElements())
    factor = Integer(card).factor()
    if list(factor)[-1][0] > 2*(n):
        print("COUNTER-EXAMPLE")
        print(deco)
        print(card)
        print(factor)
        print()
    i+=1
    if i %100 ==0:
        print(deco)
        print(card)
        print(factor)
        print()
"""
