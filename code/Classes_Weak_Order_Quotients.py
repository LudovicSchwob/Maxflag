#!/usr/bin/env python
# coding: utf-8

# # Header

# In[1]:


from sage.graphs.independent_sets import IndependentSets
from sage.plot.arc import Arc
from collections import deque
import random
import time


from IPython.core.magic import register_cell_magic

@register_cell_magic
def skip(line, cell):
    return


# # Classes

# ## SingleArc

# In[2]:


#The "SingleArc" class is an object that represents a single arc in an arc diagram 
class SingleArc:
    """
    A SingleArc is initialized with:
        an int n: in which symmetric group it belongs,
        an int i: its starting point,
        an int j: its ending point,
        a list of int L: all the values in ]i,j[ that the SingleArc passes over
        a list of int R: all the values in ]i,j[ that the SingleArc passes under
    """

    def __init__(self, i,j,L,R,n):
        self.n = n
        self.i = i
        self.j = j
        self.L = L
        self.R = R

    def __hash__(self):
        return hash((self.n,self.i,self.j,tuple(self.L)))


    #A SingleArc that starts at i, ends at j, passes over L and passes under R has a string representation of the form:
    #"(L,i,j,R)_n"

    def __str__(self):
        return f"({self.L},{self.i},{self.j},{self.R})_{self.n}"

    def __repr__(self):
        return f"({self.L},{self.i},{self.j},{self.R})_{self.n}"

    #Two SingleArcs are equal if they share the same attributes
    def __eq__(self,arc):
        #You can only compare SingleArcs with another SingleArc
        if not isinstance(arc, SingleArc):
            return NotImplemented
        return self.n == arc.n and self.i == arc.i and self.j == arc.j and self.L == arc.L and self.R == arc.R

    """
        A SingleArc is greater than another SingleArc in the weak order sense if its set of inversions contains the other's.
    """

    def __ge__(self,other):
        #print(self,self.invSet(),"vs",other,other.invSet())
        return set(other.invSet()).issubset(set(self.invSet()))




    ########################################
    #### Class Methods #####################
    ########################################

    """
    fromArc

    This class method constructs a SingleArc object from an "arc" tuple of the form "(L,i,j,R)" with "n" providing the dimension

    Return:
        SingleArc with attributes (i,j,L,R,n)
    """

    @classmethod
    def fromArc(cls,arc1,n):
        return cls(arc1[1],arc1[2],arc1[0],arc1[3],n)

    """
    isListSymmetric

    This class method determines whether the given list of SingleArc(s) is symmetric in the sense that for every SingleArc
    in the list, the list also contains its flip.

    Return:
        True -- if the list is symmetric
        False -- if not
    """

    @classmethod
    def isListSymmetric(cls,list_of_arcs_input):
        #We work on a copy of the argument
        list_of_arcs = list(list_of_arcs_input)
        #Obviously, if the input list has an odd number of SingleArc, it can't be symmetric
        if len(list_of_arcs)%2 == 1:
            return False

        k=0
        while k < len(list_of_arcs)-1:
            #print("k =",k)
            #rint(len(list_of_arcs))
            #The symm variable keeps in track whether each encountered SingleArc in the list has had it's flip present
            #The current arc hasn't seen it's flip yet
            symm = False
            curr_arc = list_of_arcs[k]

            #We iterate through the rest of the list
            for l in range(k+1,len(list_of_arcs)):
                #We check whether the current arc is the flip of another arc with the flipped method of this class
                if curr_arc == list_of_arcs[l].flipped():
                    #If we find its flip, we can remove it from the list and we can stop the search for the current SingleArc
                    symm = True
                    list_of_arcs.pop(l)
                    break
            #Here we have searched through the whole list to find the flip of the current SingleArc. If we have not found it, then the list is not symmetric
            if symm == False:
                return False
            k+=1
        #At this point we have gone through the whole list without finding any unpaired SingleArc. The list is thus symmetric
        return True


    """
        allArcs

        This class method provides a list of all the possible SingleArc(s) in an ArcDiag of size n.

        Return:
            list of SingleArc(s)
    """

    @classmethod
    def allArcs(cls,n,ideal_list = [],as_arc_list=False,essential=False):
        #Initialize the list
        list_of_arcs = []
        #We first go through every possible endpoints combinations
        if as_arc_list:
            for i in range(1,n):
                for j in range (i+1,n+1):
                    #To get every possible SingleArc with endpoints (i,j), we go over every possible configuration of over/under for the interior points.
                    #If there's no interior point, simply add the arc
                    if j-i < 2:
                        if not essential:
                            list_of_arcs.append(([],i,j,[]))
                    else:
                        #We create a list of all subsets of interior points
                        intr_points = [x for x in range(i+1,j)]
                        subsetsb = [[]]
                        for intr_point in intr_points:
                            subsetsb += [s+[intr_point] for s in subsetsb]

                        for subsetb in subsetsb:
                            list_of_arcs.append((subsetb,i,j,R))
        else:

            for i in range(1,n):
                for j in range (i+1,n+1):
                    #To get every possible SingleArc with endpoints (i,j), we go over every possible configuration of over/under for the interior points.
                    #If there's no interior point, simply add the arc
                    if j-i < 2:
                        if not essential:
                            new_single_arc = cls.fromArc(([],i,j,[]),n)
                            if new_single_arc not in ideal_list:
                                list_of_arcs.append(new_single_arc)
                    else:
                        #We create a list of all subsets of interior points
                        intr_points = [x for x in range(i+1,j)]
                        subsetsb = [[]]
                        for intr_point in intr_points:
                            subsetsb += [s+[intr_point] for s in subsetsb]

                        for subsetb in subsetsb:
                            #For every subset we add a SingleArc to the list where L contains the points in the subset and R contains the other interior points
                            R = [x for x in range(i+1,j) if x not in subsetb]
                            new_single_arc = cls.fromArc((subsetb,i,j,R),n)
                            if new_single_arc not in ideal_list:
                                list_of_arcs.append(new_single_arc)
        return list_of_arcs


    """
        forcingOrder

        This class method returns the forcing order of the weak order of size n

        Return:
            Poset object
    """



    @classmethod
    def forcingOrder(cls,n):
        data = {}
        for single_arc in cls.allArcs(n):
            if single_arc.i > 1 or single_arc.j<n:
                data[single_arc] = single_arc.forcingCovers()
        return Poset(data)



    """
        posetOfJI

        This class method returns the poset of irreducibles of the weak order of size n (i.e. the poset of every possible SingleArc in an ArcDiag of size n ordered by inversion set)

        Return:
            Poset object
    """

    @classmethod
    def posetOfJI(cls,n):
        return Poset((cls.allArcs(n), lambda x, y: x<=y))

    """
        posetOfMI

        This class method returns the poset of meet irreducibles of the weak order of size n (i.e. the poset of every possible SingleArc in an ArcDiag of size n reverse ordered by inversion set)

        Return:
            Poset object
    """

    @classmethod
    def posetOfMI(cls,n):
        #for single_arc in cls.allArcs(n):
        #    print(single_arc)
        print(cls.allArcs(n))
        return Poset((cls.allArcs(n), lambda x,y: x.meetLE(y)))








    ########################################
    #### Methods ###########################
    ########################################

    """
        forcingCovers

        This method returns the list of SingleArc(s) the cover this SingleArc in the forcing order

        Return:
            list ofSingleArc(s)
    """

    def forcingCovers(self):
        covers = []
        if self.i > 1:
            covers.append(SingleArc.fromArc(([self.i]+self.L,self.i-1,self.j,self.R),self.n))
            covers.append(SingleArc.fromArc((self.L,self.i-1,self.j,[self.i]+self.R),self.n))
        if self.j<self.n:
            covers.append(SingleArc.fromArc((self.L+[self.j],self.i,self.j+1,self.R),self.n))
            covers.append(SingleArc.fromArc((self.L,self.i,self.j+1,self.R+[self.j]),self.n))
        return covers

    """
        meetLE

        This method returns whether this meet SingleArc's associated permutation is below the input meet SingleArc's associated permutation in the weak order

        Return:
            True -- this meet SingleArc is below the input meet SingleArc
            False -- otherwise
    """

    def meetLE(self,other):
        #print("*******************")
        #print(self)
        #print(other)
        return set(other.meetInvSet()).issubset(set(self.meetInvSet()))

    """
        meetPermu

        This method returns the permutation associated to the meet arc diagram corresponding to this SingleArc

        Return:
            permutation (list of int {1,...,n})

    """

    def meetPermu(self):
        permu = []
        for val in range(self.n-self.j):
            permu.append(self.n - val)
        permu += self.L[::-1]
        permu.append(self.i)
        permu.append(self.j)
        permu += self.R[::-1]
        for val in range(1,self.i):
            permu.append(self.i-val)
        return permu



    """
    isSubarc

    This method checks if this SingleArc is a subarc of the SingleArc in the output.

    Return:
        True  -- if this SingleArc is a subarc of arc2,
        False -- if not
    """
    def isSubarc(self,arc2):

        #This checks if the interval of arc2 ([i,j]) contains the interval of this SingleArc. Returns "False" if not.
        if self.i < arc2.i or self.j > arc2.j:
            return False

        #This checks if all values this SingleArc passes over are also passed over by arc2. Returns "False" if not.

        for val in self.L:
            if val not in arc2.L:
                return False

        #This checks if all values this SingleArc passes under are also passed under by arc2

        for val in self.R:
            if val not in arc2.R:
                return False

        #At this point all checks are passed and this SingleArc is indeed a subarc of arc2
        return True


    """
    invSet

    This method returns the inversion set of the permutation P whose arc diagram is this SingleArc.

    Return:
        list of tuples of the form (j,i) if the pair i,j is an inversion of P
    """

    def invSet(self):

        inv_set = []

        #The endpoints of this SingleArc are an inversion of P
        inv_set.append((self.j,self.i))

        #All points in L are in an inversion with i and all points in R that are to its left
        for point in self.L:
            inv_set.append((point,self.i))
            for point2 in self.R:
                if point2 < point:
                    inv_set.append((point,point2))

        #All points in R are in an inversion with j
        for point in self.R:
            inv_set.append((self.j,point))


        return inv_set

    """
    meetInvSet

    This method returns the inversion set of the permutation P whose meet arc diagram is this SingleArc.

    Return:
        list of sets of the form (j,i) if the pair i,j is an inversion of P
    """

    def meetInvSet(self):
        #Initialize the set of inversions
        #We obtain the permutation that has this SingleArc as meet diagram




        return Perm(self.meetPermu()).invSet()


    """
    isCompatible

    This method returns whether this SingleArc and the input SingleArc are crossing or non-crossing.

    Return:
        True -- if the two SingleArcs are non-crossing
        False -- if the two SingleArcs cross
    """

    def isCompatible(self,arc2):

        #The two SingleArcs cross if they start at the same point and/or end at the same point
        if self.i == arc2.i or self.j == arc2.j:
            return False

        #The two SingleArcs cannot cross if one of them ends before or where the other starts
        elif self.i >= arc2.j or self.j <= arc2.i:
            return True

        else:
            #We start the comparison of the two arcs at the furthest start point and end it at the nearest end point
            start_point = max(self.i,arc2.i)
            end_point = min(self.j,arc2.j)

            #We will use the variable pos_rel (relative position) as such: 
            ## pos_rel = 1 if arc2 is over this SingleArc
            ## pos_rel = -1 if arc2 is under this SingleArc

            #In the first comparison we establish if at first arc2 is over or under this SingleArc
            if start_point == self.i:
                if start_point in arc2.L:
                    pos_rel = 1
                else:
                    pos_rel = -1
            else:
                if start_point in self.L:
                    pos_rel = -1
                else:
                    pos_rel = 1

            #Now, this SingleArc and arc2 cross only if they switch realtive position. We test it at each point
            for point in range(start_point+1,end_point):
                if (point in self.L) and (point in arc2.R):
                    if pos_rel == 1:
                        return False
                elif (point in self.R) and (point in arc2.L):
                    if pos_rel == -1:
                        return False

            #Finally, we look at the last point of comparison
            if end_point == self.j:
                if pos_rel == 1 and end_point in arc2.R:
                    return False
                elif pos_rel == -1 and end_point in arc2.L:
                    return False
            else:
                if pos_rel == 1 and end_point in self.L:
                    return False
                elif pos_rel == -1 and end_point in self.R:
                    return False

            ##If all checks are done then they are non-crossing
            return True

    """
    flipCut

    This method returns a new SingleArc that is this SingleArc where a single flip cut was applied to from the input startpoint to the input endpoint

    Return:
        SingleArc
    """

    def flipCut(self,start,end):
        #We iterate through this SingleArc's L and R sets and switch the position of the arc around points that are in the input interval
        new_L_set = []
        new_R_set = []
        for value in self.L:
            if value>start and value<end:
                new_R_set.append(value)
            else:
                new_L_set.append(value)
        for value in self.R:
            if value>start and value<end:
                new_L_set.append(value)
            else:
                new_R_set.append(value)

        new_L_set.sort()
        new_R_set.sort()

        return SingleArc(self.i,self.j,new_L_set,new_R_set,self.n)




    """
    isLinked

    This method returns whether this SingleArc is "directly linked" (Marin, Novelli, Pilaud 2025+) to the argument SingleArc

    Return:
        True -- if the two SingleArc(s) are "linked"
        False -- if not
    """



    def isLinked(self,arc2):
        #Two SingleArc(s) are "directly linked" if they share interior points.
        #We start by establishing the interval of interior points

        start_point = max(self.i,arc2.i)
        end_point = min(self.j,arc2.j)

        #They are not linked only if they don't share interior points
        if end_point - start_point < 2:
            return False

        #If they share interior points, they are linked
        return True





    """ 
    The following implementation of isLinked was following the old conjecture that is now proven to be false
    """


    def isLinked_old(self, arc2):
        #Two SingleArc(s) are "directly linked" if they share interior points and they agree on all of them (or disagree on all of them)
        #We start by establishing the interval of interior points that they share

        start_point = max(self.i,arc2.i)
        end_point = min(self.j,arc2.j)

        #If they don't share interior points, they are automatically not "directly linked"
        if end_point - start_point < 2:
            return False

        #print("Interval:", start_point,end_point)

        #Since they can either agree on every point OR disagree on every point, the first interior point will fix 
        #a value of agreeability
        if (start_point+1 in self.L and start_point+1 in arc2.L) or (start_point+1 in self.R and start_point+1 in arc2.R):
            agree = True
        else:
            agree = False

        #print("agree:",agree)

        #Now every other interior point has to follow the same agreeability
        for k in range(start_point+2,end_point):
            if not agree and ((k in self.L and k in arc2.L) or (k in self.R and k in arc2.R)):
                return False

            elif agree and ((k in self.L and k in arc2.R) or (k in self.R and k in arc2.L)):
                return False

        #At this point all interior points have been verified and they all followed the same agreeability.
        #The two SingleArc(s) are thus linked
        return True









    """
    flipped

    This method returns a flipped copy of this SingleArc

    Return:
        SingleArc with attributes (i,j,R,L,n)
    """

    def flipped(self):
        return SingleArc(self.i,self.j,self.R,self.L,n)



    ##
    #Graph Methods
    ##

    """
    sepArcs

    This method returns a list of the path this SingleArc follows separated into semicricles

    Return:
        list of tuples of the form {[(s,e,h), ...]} where
            s is an int at the start point of the semi circle
            e is an int at the end of the semi circle
            h is an int taking values {0,1,2}. 
                0 when the semicircle is over the abscissa
                1 when the semicircle is under the abscissa
                2 when the "semicircle" lies on the abscissa
    """

    def sepArcs(self):

        #When this SingleArc is has length 1 (i.e. j = i+1), then it lies on the abscissa
        if self.j - self.i == 1:
            return [(self.i,self.j,2)]

        else:
            ##We initialize the list of semicircles
            list_demi_arcs = []

            #The starting point of the first semicircle is i
            s_point = self.i

            #Now we determine whether this semicircle is over or under the abscissa
            if len(self.L)>0 and self.L[0] == self.i+1:
                #In this case, i+1 is in L meaning this SingleArc passes over i+1
                curr_half = 0
            else:
                #In this case, i+1 is not in L meaning it's in R and this SingleArc passes under i+1
                curr_half = 1

            #Now we iterate through all the points in the interval of this SingleArc.
            #Whenever points l and l+1 don't belong in the same set L or R, this SingleArc crosses the abscissa at
            #the midpoint (2l+1)/2.
            for l in range(self.i+2,self.j):
                if (curr_half == 0 and l in self.R) or (curr_half == 1 and l in self.L):
                    #Here, l does not belong in the same set as the previous point

                    #The end point of the current semicircle is the midpoint between l and l-1
                    e_point = l-.5

                    #We add this completed semicircle to the list
                    list_demi_arcs.append((s_point,e_point,curr_half))

                    #The new semicircle starts where the current one ends
                    s_point= e_point

                    #This new semicircle will be on the opposite half of the one just completed
                    curr_half = 1-curr_half

            #We have reached j (the end of this SingleArc). The current (and last) semicircle thus ends at j.
            list_demi_arcs.append((s_point,self.j,curr_half))

            return list_demi_arcs

    ##
    #Graph Methods
    ##

    """
    graph

    This method produces a Graph object that is a drawing of this SingleArc

    return:
        a sage Graph object representing this SingleArc
    """

    def graph(self,shift=0,anchor = 0, direction =0):
        #We start with the points. One point for each value {1,...,n}
        point_list = []
        n = self.n
        for l in range(1,self.n+1):
            point_list.append([l+shift,0])

        #We initialize the graph object
        graph = point(point_list,color = 'black')

        #If there is an anchor to display, put it in green
        if anchor >0 and anchor <= self.n:
            graph+= point([anchor+shift,0])

        if direction == -1:
            #Display an arrow to the left
            graph+=arrow((anchor+shift+.5,-1),(anchor+shift-.5,-1))

        if direction == 1:
            #Display an arrow to the right
            graph+=arrow((anchor+shift-.5,-1),(anchor+shift+.5,-1))

        #Now we add the arc each semicircle at a time
        semi_circles = self.sepArcs()
        if semi_circles[0][2] == 2:
            #In this case, this SingleArc is a line between two adjacent points
            graph += line([(semi_circles[0][0]+shift,0),(semi_circles[0][1]+shift,0)], color = 'red',thickness = 3)
        else:
            for semi_circle in semi_circles:
                curr_half = semi_circle[2]
                #We use the sage arc graph object for the semicircles
                graph += arc(((semi_circle[0]+semi_circle[1])/2+shift,0),(semi_circle[1]-semi_circle[0])/2, sector = (0+curr_half*pi,pi+curr_half*pi),thickness=3, color = 'red')

        #Finally, we draw a box around the drawing
        graph += line([(0+shift,n/2+1),(n+1+shift,n/2+1)])
        graph += line([(0+shift,-(n/2+1)),(n+1+shift,-(n/2+1))])
        graph += line([(0+shift,n/2+1),(0+shift,-(n/2+1))])
        graph += line([(n+1+shift,n/2+1),(n+1+shift,-(n/2+1))])
        graph.axes(False)
        return graph




# ## Perm

# In[4]:


#The "Perm" class is an object that represents a permutation. It will be useful for developping methods
class Perm:
    """
    A Perm is initialized with:
        a list perm: the permutation on  {1,...,n}
    """

    def __init__(self,perm):
        self.perm = perm

    """
    Hash function
    """

    def __hash__(self):
        return hash(tuple(self.perm))

    """
    Two Perm is greater than an other (in terms of the weak order) if its set of inversions contains the set of inversion of the other Perm
    """

    def __ge__(self,other):
        return set(other.invSet()).issubset(set(self.invSet()))

    def __gt__(self,other):
        selfset = set(self.invSet())
        otherset = set(other.invSet())
        return otherset.issubset(selfset) and otherset != selfset

    def __str__(self):
        return f'{list(self.perm)}'

    def __repr__(self):
        return f'{list(self.perm)}'

    #Two Perm(s) are equal if they share the same permutation
    def __eq__(self,other):
        #You can only compare Perm(s) with another Perm
        if not isinstance(other, Perm):
            return NotImplemented

        return self.perm == other.perm





    ########################################
    #### Class Methods #####################
    ########################################


    """
    invSet

    This method gives the set of inversions of this Perm as a list

    Return:
        list of sets of the form {i,j} if the pair i,j is an inversion of this Perm
    """

    def invSet(self):
        #We will obtain the inversions by looking at the permutation one position at a time, and for each position, we look at a list of the values before the current position. Each time the current value is smaller than a value in
        #we add this inversion to the list.
        inv_set = []
        visited_vals = []
        for curr_val in self.perm:
            for visited_val in visited_vals:
                if curr_val<visited_val:
                    inv_set.append((visited_val,curr_val))

            visited_vals.append(curr_val)
        return inv_set

    """
    nonCrossingArcDiag

    This method gives the non-crossing ArcDiag corresponding to this Perm

    Return:
        ArcDiag of dimension n corresponding to this Perm
    """

    def nonCrossingArcDiag(self):
        #We first obtain the dimension of the ArcDiag
        n = len(self.perm)

        #We now construct the list of arcs as a list of tuples
        list_arc = []

        #Each arc in an ArcDiag corresponding to a Perm comes from a descent from that Perm
        for i in range(n-1):
            if self.perm[i]>self.perm[i+1]:
                #Here we are at a descent of this Perm. We initialize the arc for this descent
                curr_arc = ([],self.perm[i+1],self.perm[i],[])

                #If the descent is between two adjacent values, this arc is complete
                if self.perm[i] - self.perm[i+1] == 1:
                    #print("reached")
                    list_arc.append(curr_arc)
                    continue

                #We read the shortest sliced list between self.perm[:i] and self.perm[i+1:]
                if i < int((n+1)/2)-1:
                    #print((self.perm[i],self.perm[i+1]),"read left")
                    for value in self.perm[:i]:
                        if value > self.perm[i+1] and value < self.perm[i]:
                            #print("value:",value)
                            #Here value is between the two values of the descent of the arc and before the descent
                            #in the permutation. It is sent to L
                            curr_arc[0].append(value)
                    #L now contains all its values and is sorted
                    curr_arc[0].sort()
                    #R will now contain all values between the two values of the descent that are not in L

                    #There are two cases to consider. Either L is empty and we put all values between the two values
                    #in the descent in R, or L is not empty and we only put the missing values in R.
                    if len(curr_arc[0]) == 0:
                        #Here L is empty so every intermediate value are sent to R
                        for value in range(self.perm[i+1]+1,self.perm[i]):
                            curr_arc[3].append(value)
                    else:
                        #Here L is not empty. We go through every intermediate value of the arc
                        for value in range(self.perm[i+1]+1,self.perm[i]):
                            #We declare a variable to_add which keeps track whether the current value needs to be added or not. It begins as True
                            to_add = True
                            #Now go through L which is sorted
                            for valueL in curr_arc[0]:
                                #Going through L in order, if we meet the value, then we don't add it in R and we proceed to the next value.
                                #If we meet a value greater than the current value, this means it's not present in L and needs to be added in R.
                                #If we go through L and have not seen the current value or a value greater, in L, then we also need to add it to R
                                if valueL>value:
                                    break
                                elif valueL == value:
                                    to_add = False
                                    break
                                else:
                                    continue
                            if to_add:
                                curr_arc[3].append(value)
                else:
                    #print((self.perm[i],self.perm[i+1]),"read right")
                    for value in self.perm[i+1:]:
                        if value > self.perm[i+1] and value < self.perm[i]:
                            #print("value:",value)
                            #Here value is between the two values of the descent of the arc and before the descent 
                            #in the permutation. It is sent to R
                            curr_arc[3].append(value)
                    #R now contains all its values and is sorted
                    curr_arc[3].sort()
                    #L will now contain all values between the two values of the descent that are not in R

                    #There are two cases to consider. Either R is empty and we put all values between the two values
                    #in the descent in L, or R is not empty and we only put the missing values in L.
                    if len(curr_arc[3]) == 0:
                        #Here R is empty so every intermediate value are sent to L
                        for value in range(self.perm[i+1]+1,self.perm[i]):
                            curr_arc[0].append(value)
                    else:
                        #Here L is not empty. We go through every intermediate value of the arc
                        for value in range(self.perm[i+1]+1,self.perm[i]):
                            #We declare a variable to_add which keeps track whether the current value needs to be added or not. It begins as True
                            to_add = True
                            #Now go through L which is sorted
                            for valueR in curr_arc[3]:
                                #Going through L in order, if we meet the value, then we don't add it in R and we proceed to the next value.
                                #If we meet a value greater than the current value, this means it's not present in L and needs to be added in R.
                                #If we go through L and have not seen the current value or a value greater, in L, then we also need to add it to R
                                if valueR>value:
                                    curr_arc[0].append(value)
                                    to_add = False
                                    break
                                elif valueR == value:
                                    to_add = False
                                    break
                                else:
                                    continue
                            if to_add:
                                curr_arc[0].append(value)

                #The arc is complete and added to the list
                list_arc.append(curr_arc)
        #All descents of the permuatation have had their arcs added to the list. The list is complete.
        return ArcDiag.fromArcList(list_arc,n)

    """
        meetNonCrossing

        This method gives the blue non-crossing ArcDiag corresponding to this Perm

        Return:
            ArcDiag of dimension n corresponding to this Perm's meet irreducible representation
    """

    def meetNonCrossing(self):
        #We first obtain the dimension of the ArcDiag
        n = len(self.perm)

        #We now construct the list of arcs as a list of tuples
        list_arc = []

        #Each arc in an ArcDiag corresponding to a Perm comes from an ascent from that Perm
        for i in range(n-1):
            if self.perm[i]<self.perm[i+1]:
                #Here we are at an ascent of this Perm. We initialize the arc for this ascent
                curr_arc = ([],self.perm[i],self.perm[i+1],[])

                #If the ascent is between two adjacent values, this arc is complete
                if self.perm[i+1] - self.perm[i] == 1:
                    #print("reached")
                    list_arc.append(curr_arc)
                    continue

                #We read the shortest sliced list between self.perm[:i] and self.perm[i+1:]
                if i < int((n+1)/2)-1:
                    #print((self.perm[i],self.perm[i+1]),"read left")
                    for value in self.perm[:i]:
                        if value > self.perm[i] and value < self.perm[i+1]:
                            #print("value:",value)
                            #Here value is between the two values of the descent of the arc and before the descent
                            #in the permutation. It is sent to L
                            curr_arc[0].append(value)
                    #L now contains all its values and is sorted
                    curr_arc[0].sort()
                    #R will now contain all values between the two values of the descent that are not in L

                    #There are two cases to consider. Either L is empty and we put all values between the two values
                    #in the descent in R, or L is not empty and we only put the missing values in R.
                    if len(curr_arc[0]) == 0:
                        #Here L is empty so every intermediate value are sent to R
                        for value in range(self.perm[i]+1,self.perm[i+1]):
                            curr_arc[3].append(value)
                    else:
                        #Here L is not empty. We go through every intermediate value of the arc
                        for value in range(self.perm[i]+1,self.perm[i+1]):
                            #We declare a variable to_add which keeps track whether the current value needs to be added or not. It begins as True
                            to_add = True
                            #Now go through L which is sorted
                            for valueL in curr_arc[0]:
                                #Going through L in order, if we meet the value, then we don't add it in R and we proceed to the next value.
                                #If we meet a value greater than the current value, this means it's not present in L and needs to be added in R.
                                #If we go through L and have not seen the current value or a value greater, in L, then we also need to add it to R
                                if valueL>value:
                                    break
                                elif valueL == value:
                                    to_add = False
                                    break
                                else:
                                    continue
                            if to_add:
                                curr_arc[3].append(value)
                else:
                    #print((self.perm[i],self.perm[i+1]),"read right")
                    for value in self.perm[i+1:]:
                        if value > self.perm[i] and value < self.perm[i+1]:
                            #print("value:",value)
                            #Here value is between the two values of the ascent of the arc and after the ascent 
                            #in the permutation. It is sent to R
                            curr_arc[3].append(value)
                    #R now contains all its values and is sorted
                    curr_arc[3].sort()
                    #L will now contain all values between the two values of the descent that are not in R

                    #There are two cases to consider. Either R is empty and we put all values between the two values
                    #in the descent in L, or R is not empty and we only put the missing values in L.
                    if len(curr_arc[3]) == 0:
                        #Here R is empty so every intermediate value are sent to L
                        for value in range(self.perm[i]+1,self.perm[i+1]):
                            curr_arc[0].append(value)
                    else:
                        #Here L is not empty. We go through every intermediate value of the arc
                        for value in range(self.perm[i]+1,self.perm[i+1]):
                            #We declare a variable to_add which keeps track whether the current value needs to be added or not. It begins as True
                            to_add = True
                            #Now go through L which is sorted
                            for valueR in curr_arc[3]:
                                #Going through L in order, if we meet the value, then we don't add it in R and we proceed to the next value.
                                #If we meet a value greater than the current value, this means it's not present in L and needs to be added in R.
                                #If we go through L and have not seen the current value or a value greater, in L, then we also need to add it to R
                                if valueR>value:
                                    curr_arc[0].append(value)
                                    to_add = False
                                    break
                                elif valueR == value:
                                    to_add = False
                                    break
                                else:
                                    continue
                            if to_add:
                                curr_arc[0].append(value)

                #The arc is complete and added to the list
                list_arc.append(curr_arc)
        #All ascents of the permutation have had their arcs added to the list. The list is complete.
        return ArcDiag.fromArcList(list_arc,n)

    ##
    #Graph Methods
    ##


    """
    graph

    This method produces a Graph object that is a drawing of this Perm's non crossing arc diagram

    return:
        a sage Graph object representing this Perm's non crossing arc diagram
    """

    def graph(self,shift=0,anchor = 0,direction = 0):
        arc_diag = self.nonCrossingArcDiag()
        return arc_diag.graph(shift=shift,anchor = anchor,direction=direction)

    """
    meetGraph

    This method produces a Graph object that is a drawing of this Perm's non crossing arc diagram

    return:
        a sage Graph object representing this Perm's non crossing arc diagram
    """

    def meetGraph(self,shift=0):
        meet_arc_diag = self.meetNonCrossing()
        return meet_arc_diag.graph(shift,arc_color = 'blue')






# ## ArcDiag

# In[5]:


  #The "ArcDiag" class is an object that represents a (not necessarily non-crossing) arc diagram
class ArcDiag:
    """
    An ArcDiag is initialized with:
        an int n: it's dimension
        a list of SingleArc(s) arc_list: the arcs of this ArcDiag
    """

    def __init__(self,arc_list,n):
        self.n = n
        self.arc_list = arc_list

    #An ArcDiag of dimension n and k arcs in arc_list has a string representation of the form:
    #"An arc diagram of dimension {n} and arcs: {arc_list[0]}, {arc_list[1]}, ..., {arc_list[k-1]}"
    def __str__(self):
        return f'An arc diagram of dimension {self.n} and arcs: {" ".join(f"{self.arc_list[i]}" for i in range(len(self.arc_list)))}'

    #Two ArcDiag(s) are the same if they have the same SingleArc(s) and the same n
    def __eq__(self,other):
        if self.n != other.n:
            return False
        other_copy = [x for x in other.arc_list]
        for x in self.arc_list:
            if x not in other_copy:
                return False
            other_copy.remove(x)
        if len(other_copy)>0:
            return False
        return True

    ########################################
    #### Class Methods #####################
    ########################################

    """
    fromArcList

    This method constructs an ArcDiag object from a list of "arc" tuple of the form "[(L,i,j,R),...]" with "n" providing the dimension

    Return:
        ArcDiag of dimension {n} and arcs {[(L,i,j,R),...]}
    """

    @classmethod
    def fromArcList(cls,arcs,n):
        #We construct the arc_list by transforming each "arc" tuple into a SingleArc with the fromArc class method 
        arc_list = []
        for i in range(len(arcs)):
            arc_list.append(SingleArc.fromArc(arcs[i],n))
        return cls(arc_list,n)

    ########################################
    #### Methods ###########################
    ########################################

    """
    listOfArcs

    This method returns the list of all arcs of this ArcDiag in the tuple representation

    Return:
        A list of tuples (each tuple is an arc of this ArcDiag)
    """

    def listOfArcs(self):
        #In this ArcDiag's arc_list, all the arcs are SingleArc objects. Now we want them in the form of tuples
        list_of_arcs = []

        for single_arc in self.arc_list:
            list_of_arcs.append((single_arc.L,single_arc.i,single_arc.j,single_arc.R))

        list_of_arcs.sort(key = lambda x : x[1])

        return list_of_arcs


    """
    isNonCrossing

    This method returns whether this ArcDiag is non-crossing or not

    Return:
        True -- if this ArcDiag is non-crossing
        False -- if this ArcDiag is not non-crossing
    """

    def isNonCrossing(self):
        #An ArcDiag is non-crossing if and only if its SingleArc(s) are pairwise compatible [Reading, arXiv:1405.6904]
        #We will test the compatibility of each arc with every other arc

        for i in range(len(self.arc_list)):
            #Compatibility is symmetric
            for j in range(i+1,len(self.arc_list)):

                #To test the compatibility of two SingleArc(s) we use the isCompatible method of the SingleArc class
                if not self.arc_list[i].isCompatible(self.arc_list[j]):
                    #Once we encounter a pair that is not compatible, we know the ArcDiag is not non-crossing
                    return False

        #At this point we have tested for each pair of SingleArc(s) and we can conclude the ArcDiag is non-crossing
        return True

    """
    sortedArcListStart

    This method sorts the arc list of this ArcDiag by their start points.

    Void method
    """

    def sortedArcListStart(self):
        self.arc_list.sort(key = lambda x : x.i)


    """
    components

    This method returns the components [Reading, arXiv:1405.6904] of this ArcDiag (granted it is non-crossing). A component is a couple of lists. The first list consists of all the endpoints of the SingleArc(s) in the component. The second list consists of all points this compoenent goes "under"


    Return:
        a list of components
    """

    def components(self):
        #Components are defined only on non-crossing ArcDiag(s).
        if not self.isNonCrossing():
            raise Exception("Not non-crossing")    

        #A component of this ArcDiag is a maximal list of points (p_1 < ... < p_k) where there exists an arc starting at p_i ending at p_(i+1) for all i in {1,...,k-1}.

        #The list "points" will keep track of the points which are not part of any component and will be components on their own
        points = [k for k in range(1,self.n+1)]

        #We first sort the SingleArc(s) of this ArcDiag by their start points
        self.sortedArcListStart()

        #We initialize the list of components
        comp_list = []

        #We iterate through all the arcs
        for single_arc in self.arc_list:
            #We need to track whether that SingleArc's endpoints were added to an already existing component
            added = False
            #We first look whether that SingleArc starts where an earlier SingleArc of this ArcDiag ended
            for comp in comp_list:
                if single_arc.i == comp[0][-1]:
                    #If it is the case, we add its end point to the component
                    comp[0].append(single_arc.j)
                    #We also add all the values that SingleArc "goes under" to the second list of the component
                    for k in range(len(single_arc.L)):
                        comp[1].append(single_arc.L[k])
                    #We also remove its end point from the list of points with no component
                    points.remove(single_arc.j)
                    #This SingleArc's endpoints were added to an already existing component
                    added = True
                    break

            #In the case where the endpoints of this SingleArc where not added to an already existing component, we create a new component
            if not added:
                comp_list.append(([single_arc.i,single_arc.j],single_arc.L))
                #Both endpoints were added to a component so they won't each be in a component on their own
                points.remove(single_arc.i)
                points.remove(single_arc.j)

        #Once we are done looking at all the SingleArc(s), the points that were not added yet are all each on a component of their own
        for point in points:
            comp_list.append(([point],[]))

        comp_list.sort(key = lambda x : x[0])
        return comp_list




    """
    perm

    This method returns the permutation associated to this ArcDiag (granted it is non-crossing)

    Return:
        permutation on {1,...,n} (NOT a perm object)
    """
    def perm(self):
        #First, if this ArcDiag has no SingleArc then the corresponding permutation is the identity
        if len(self.arc_list) == 0:
            return [k for k in range(1,self.n+1)]


        #We start by obtaining all components [Reading, arXiv:1405.6904] of this ArcDiag
        comp_list = self.components()
        #print(comp_list)

        #We initialize a list that will ultimately contain all the components in the correct order from lowest to highest.
        order_list = []

        #We iterate through all components of this ArcDiag
        for comp in comp_list:
            #To find its position in the order, we iterate through the ordered list
            i = 0
            for comp_2 in order_list:
                if comp[0][0] in comp_2[1]:
                    #The first component in the ordered list that goes over the current component stops the iterating as we have found the correct index
                    break
                i+=1

            order_list.insert(i,comp)

        #We now initiate the permutation
        perm = []
        #print(order_list)

        #To obtain the correct permutation, we read one by one each component in reverse order, starting from the first component in the ordered list
        for comp in order_list:
            comp[0].reverse()
            perm+=comp[0]

        return perm

    """
    meetPerm

    This method returns the permutation associated to this meet ArcDiag (granted it is non-crossing)

    Return:
        permutation on {1,...,n} (NOT a perm object)
    """
    def meetPerm(self):
        #First, if this ArcDiag has no SingleArc then the corresponding permutation is the inverse identity
        if len(self.arc_list) == 0:
            return [self.n+1-k for k in range(1,self.n+1)]


        #We start by obtaining all components [Reading, arXiv:1405.6904] of this ArcDiag
        comp_list = self.components()
        comp_list.sort(key = lambda x : x[0][-1])
        comp_list.reverse()
        #print(comp_list)

        #We initialize a list that will ultimately contain all the components in the correct order from lowest to highest.
        order_list = []

        #We iterate through all components of this ArcDiag
        for comp in comp_list:
            #To find its position in the order, we iterate through the ordered list
            i = 0
            for comp_2 in order_list:
                if comp[0][0] in comp_2[1]:
                    #The first component in the ordered list that goes over the current component stops the iterating as we have found the correct index
                    break
                i+=1

            order_list.insert(i,comp)
        #print(order_list)

        #We now initiate the permutation
        perm = []

        #To obtain the correct permutation, we read one by one each component in reverse order, starting from the last component in the ordered list
        for comp in order_list:
            perm+=comp[0]

        return perm








    """
    isIdeal

    This method returns whether this ArcDiag is an ideal on the forcing order of S_n of not

    Return:
        True -- if this ArcDiag is an ideal on the forcing order of S_n
        False -- if not
    """

    def isIdeal(self):
        #An ArcDiag is an ideal on the forcing order of S_n if and only if its SingleArc(s) are pairwise not subarc of one another
        #We will test for each pair if one is subarc of the other

        for i in range(len(self.arc_list)):
            #The subarc relation is not symmetric so for each pair we will need to test each direction
            for j in range(i+1,len(self.arc_list)):

                #To test whether a SingleArc is subarc of another SingleArc we use the isSUbarc method of the SingleArc class
                if self.arc_list[i].isSubarc(self.arc_list[j]) or self.arc_list[j].isSubarc(self.arc_list[i]):
                    #Once we encounter a SingleArc that is subarc of another SingleArc we know this ArcDiag cannot be ideal
                    return False

        #At this point we have tested for each pair of SingleArc(s) and we can conclude the ArcDiag is ideal on the forcing order of S_n
        return True


    """
    flipped

    This method returns a flipped copy of this ArcDiag

    Return:
        ArcDiag (flipped version of this ArcDiag)
    """

    def flipped(self):
        new_arc_list = []
        for single_arc in self.arc_list:
            new_arc_list.append(single_arc.flipped())
        return ArcDiag(new_arc_list,self.n)


    """
    idealCompatible

    This method returns whether this ArcDiag has no subarcs among the SingleArc(s) of the input ArcDiag

    Return:
        False -- If at least one of the SingleArc(s) of this ArcDiag has a subarc in the input ArcDiag
        True -- Otherwise

    """
    def idealCompatible(self,ideal):
        #We simply verify for each arc of this ArcDiag whether one of the arcs of the input ArcDiag is its subarc
        for single_arc in self.arc_list:
            for ideal_arc in ideal.arc_list:
                if ideal_arc.isSubarc(single_arc):
                    return False
        return True

    ###
    # Bijection Methods #
    ###



    """
    untangle

    This method return an ArcDiag that is this ArcDiag where all crossings are replaced by "pairs of elbows"

    Return:
        ArcDiag -- That ArcDiag will be non-crossing
    """

    def untangle(self):
        #We will iterate left through right, point by point through the diagram. 
        arc_queue = []
        #As we will have to add the arcs of this ArcDiag one by one, it is a good idea to store all of its arcs in a sorted copy of it's arclist
        #They are sorted by their starting points in reverse order (We want to remove the first ones first
        arc_list_copy = sorted(self.arc_list, key = lambda x : x.i,reverse = True)

        #We will also need to keep track of the orders of the SingleArc(s) among themselves. However, the Single arcs will change shape in the process so we order the starting points instead
        order_vector = []

        #Finally, we also initialize a list of all the arcs in the untangled ArcDiag
        final_arc_list = []

        for curr_pos in range(1,self.n+1):

            #We begin by looking through all accumulated SingleArc(s) at this point and see if they crossed tho we only do this next step if there are at least two SingleArc(s) in the queue

            if len(arc_queue)>1:
                i = 0
                #We use a while loop here because the list will keep changing
                while i < len(arc_queue):

                    #We keep the start points of the first comparison SingleArc
                    start_i = arc_queue[i].i
                    end_i = arc_queue[i].j
                    #This variable keeps track of whether a change was made to the current SingleArc. If yes, since its new version will be sent at the end of the list we will not advance in the list
                    change_occured = 0

                    #Now we check where the fisrt compraisson SingleArc is placed with regard to the current point
                    if curr_pos == arc_queue[i].j:
                        prev_place = 1
                    elif curr_pos in arc_queue[i].R:
                        prev_place = 0

                    else:
                        prev_place = 2

                    #Now we look through the remaining SingleArc(s) of the list to see if any cross
                    for j in range(i,len(arc_queue)):

                        #We check where the second comparisson SingleArc is placed with regard to the current point
                        start_j = arc_queue[j].i
                        end_j = arc_queue[j].j

                        if curr_pos == arc_queue[j].j:
                            curr_place = 1
                        elif curr_pos in arc_queue[j].R:
                            curr_place = 0

                        else:
                            curr_place = 2

                        #We check whether the two SingleArc(s) cross

                        if (order_vector.index(start_i) < order_vector.index(start_j)) and curr_place < prev_place:
                            #If they do, we uncross them. We also know that they cross precisely between the last point and the current point so to uncross them we chop both SingleArc(s) into before and after
                            #the current point. We glue the lowest first half with the lowest last half and highest first half with the highest last half
                            new_low_L = []
                            new_low_R = []
                            new_high_L = []
                            new_high_R = []

                            for k in range(start_i+1,curr_pos):
                                if k in arc_queue[i].L:
                                    new_low_L.append(k)
                                else:
                                    new_low_R.append(k)

                            for k in range(start_j+1,curr_pos):
                                if k in arc_queue[j].L:
                                    new_high_L.append(k)
                                else:
                                    new_high_R.append(k)

                            for k in range(curr_pos,end_i):
                                if k in arc_queue[i].L:
                                    new_high_L.append(k)
                                else:
                                    new_high_R.append(k)

                            for k in range(curr_pos,end_j):
                                if k in arc_queue[j].L:
                                    new_low_L.append(k)
                                else:
                                    new_low_R.append(k)

                            new_low_R.sort()
                            new_low_L.sort()
                            new_high_L.sort()
                            new_high_R.sort()

                            #We remove the current SingleArc(s) and replace them with the new glued SingleArc(s) appended at the end of the list. 

                            arc_queue.pop(j)
                            arc_queue.pop(i)
                            arc_queue.append(SingleArc.fromArc((new_low_L,start_i,end_j,new_low_R),self.n))
                            arc_queue.append(SingleArc.fromArc((new_high_L,start_j,end_i,new_high_R),self.n))






                            change_occured = 1
                            break


                        elif (order_vector.index(arc_queue[j].i) < order_vector.index(arc_queue[i].i)) and prev_place < curr_place:

                            new_low_L = []
                            new_low_R = []
                            new_high_L = []
                            new_high_R = []

                            for k in range(start_i+1,curr_pos):
                                if k in arc_queue[i].L:
                                    new_high_L.append(k)
                                else:
                                    new_high_R.append(k)

                            for k in range(start_j+1,curr_pos):
                                if k in arc_queue[j].L:
                                    new_low_L.append(k)
                                else:
                                    new_low_R.append(k)

                            for k in range(curr_pos,end_i):
                                if k in arc_queue[i].L:
                                    new_low_L.append(k)
                                else:
                                    new_low_R.append(k)

                            for k in range(curr_pos,end_j):
                                if k in arc_queue[j].L:
                                    new_high_L.append(k)
                                else:
                                    new_high_R.append(k)

                            new_low_R.sort()
                            new_low_L.sort()
                            new_high_L.sort()
                            new_high_R.sort()

                            arc_queue.pop(j)
                            arc_queue.pop(i)
                            arc_queue.append(SingleArc.fromArc((new_low_L,start_j,end_i,new_low_R),self.n))
                            arc_queue.append(SingleArc.fromArc((new_high_L,start_i,end_j,new_high_R),self.n))


                            change_occured = 1
                            break

                    if change_occured== 0:
                        i+=1











            #We add arcs that start at the current point
            while len(arc_list_copy) != 0:
                #We obtain the start point of the last SingleArc. It should be weakly minimal among the start points of the remaining SingleArc(s) in the list 
                next_start_point = arc_list_copy[-1].i

                #If that start point is at the current position, we add it to the end of the queue
                if next_start_point == curr_pos:
                    arc_queue.append(arc_list_copy.pop())
                    #We also need to add it to the order vector in the correct place
                    correct_index = 0
                    for single_arc in arc_queue:
                        #Since the SingleArc we just added started at the current position, any arc lower than it has to pass under the current position point
                        if curr_pos in single_arc.R:
                            correct_index +=1

                    #At this point, we have reached the the arcs that are "higher" than the SingleArc to add
                    order_vector.insert(correct_index,arc_queue[-1].i)

                #If that start point is further, then we can stop looking as all other SingleArc(s) start weakly later than that start point
                elif next_start_point > curr_pos:
                    break

                #If that start point is earlier, then we missed including its SingleArc in an earlier interation. We raise an error
                else:
                    raise Exception("Missed a SingeArc")


            """
            print("#############")
            print(curr_pos)
            print("Order Vector:",order_vector)
            print()
            """

            #Finally, we remove the SingleArc(s) of the queue that end at the current point, and add them to the list of final SingleArc(s)
            i = 0
            while i < len(arc_queue):

                #If a SingleArc ends here, add it to the final list and remove it from the queue and from the order vector
                if arc_queue[i].j == curr_pos:
                    final_arc_list.append(arc_queue[i])
                    order_vector.remove(arc_queue[i].i)
                    arc_queue.remove(arc_queue[i])
                else:
                    i+=1
        return ArcDiag(final_arc_list,self.n)


    """
    cutAndFlip

    This method takes an ideal ArcDiag as input and transforms this ArcDiag into the output ArcDiag such that no SingleArc in the output ArcDiag has as subarc a SingleArc in the input ArcDiag.
    This is done by comparing each SingleArc of this ArcDiag with each SingleArc of the input ArcDiag. Whenever an input SingleArc A is subarc of a SingleArc B of this ArcDiag, we "cut" A out of B, flip it upside down
    and glue it back to B to obtain B'.

    Return:
        ArcDiag (possibly crossing) that has no subarc in the input ArcDiag whenever this ArcDiag is element to at least one WOQuotient in the sqame flip class as the input ArcDiag. 
    """



    def cutAndFlip(self,new_ideal):
        #We initialize the new list of arcs
        arc_list_queue = []
        new_arc_list = []
        #We deep copy each each arc of this ArcDiag's arc list into the new_arc_list
        for single_arc in self.arc_list:
            arc_list_queue.append(SingleArc(single_arc.i,single_arc.j,single_arc.L,single_arc.R,self.n))


        #We iterate through the queue of SingleArc until every SingleArc is no longer subarc
        #[We keep track of the number of iterations and of the different SingleArcs obtain in the case we ever get to a cycle (This might be impossible)]
        nb_iteration = 0
        limit = len(arc_list_queue)*len(new_ideal.arc_list)
        visited_arcs = []
        while len(arc_list_queue)>0:
            curr_single_arc = arc_list_queue.pop()
            #curr_single_arc.graph().show()
            #[We look for a possible infinite cycle]
            if (nb_iteration > limit and curr_single_arc in visited_arcs):
                """
                print("########################################")
                for vis_arc in visited_arcs:
                    print("~~~~~~~~~~~~~~~~~~~~~~~")
                    vis_arc.graph().show()
                """
                print("########################################")
                new_ideal.graph().show()
                curr_single_arc.graph().show()
                raise Exception("Possible cycle")
            #We will proceed with this iteration so we increase the number
            nb_iteration+=1
            #This variable stores whether the current SingleArc was modified and added back at the end
            was_mod=0
            #We add the SingleArc we are currently visiting into the visted list
            visited_arcs.append(SingleArc(curr_single_arc.i,curr_single_arc.j,list(curr_single_arc.L),list(curr_single_arc.R),self.n))


            for ideal_arc in new_ideal.arc_list:
                #We verify whether the current new ideal SingleArc is subarc of the current SingleArc of this ArcDiag
                if ideal_arc.isSubarc(curr_single_arc):
                    """
                    print("Change initiated:")
                    curr_single_arc.graph().show(figsize = 2)
                    ideal_arc.graph().show(figsize=2)
                    """
                    for k in range(ideal_arc.i+1,ideal_arc.j):

                        #We iterate through the internal points of the ideal SingleArc and flip the position of this ArcDiag's SingleArc with regard to the current internal point
                        if k in curr_single_arc.L:
                            curr_single_arc.R.append(k)
                            curr_single_arc.L.remove(k)
                        elif k in curr_single_arc.R:
                            curr_single_arc.L.append(k)
                            curr_single_arc.R.remove(k)
                    #We sort back the L and R lists
                    curr_single_arc.L.sort()
                    curr_single_arc.R.sort()
                    #We add the modified SingleArc back to the queue
                    arc_list_queue.append(curr_single_arc)
                    #The arc was modified and added back
                    was_mod=1
                    #We break out of the loop
                    break
            #If we reached this step without modifying the SingleArc, then this SingleArc has no subarc among the SingleArc(s) in the input ArcDiag
            if was_mod == 0:
                new_arc_list.append(curr_single_arc)
        #Here the queue was emptied so we have all the arcs
        """
        for vis_arc in visited_arcs:
            print("~~~~~~~~~~~~~~~~~~~~~~~")
            vis_arc.graph().show()
        print("########################################")
        """
        return ArcDiag(new_arc_list,self.n)




    """
    bijection (THIS DOES NOT FINISH) [([],1,3,[2]), ([3],2,5,[4]), ([],2,5,[3, 4]), ([4],3,5,[])] with ([3,4],2,5,[])

    This method takes as input an ideal ArcDiag in the same flip class of an ArcDiag whose corresponding WOQuotient has this ArcDiag as minimal element. It returns the ArcDiag of the element of the WOQuotient of the input ArcDiag corresponding to
    this ArcDiag with regard to a bijection [Marin 2025+].

    Return:
        ArcDiag that is a minimal element of the WOQuotient with ideal as the input ArcDiag
    """

    def bijection(self,new_ideal):
        new_arc_diag = self.cutAndFlip(new_ideal).untangle()
        stop = False
        while not stop:
            if not new_arc_diag.idealCompatible(new_ideal):
                new_arc_diag = new_arc_diag.cutAndFlip(new_ideal).untangle()
            else:
                stop = True
        return new_arc_diag


    """
    newBijection (THIS DOES NOT FINISH) [([],1,3,[2]) ([3],2,5,[4]) ([3, 4],2,5,[]) ([],2,5,[3, 4])] with [3, 5, 2, 4, 1]

    This method takes as input and ideal ArcDiag in the same flip class of an ArcDiag whose corresponding WOQuotient has this ArcDiag as minimal element. It returns the ArcDiag of the element of the WOQuotient of the input ArcDiag corresponding to
    this ArcDiag with regard to a bijection [Marin 2025+].

    Return:
        ArcDiag that is a minimal element of the WOQuotient with ideal as the input ArcDiag
    """

    def newBijection(self,new_ideal):
        #If this ArcDiag is already compatible with the input ArcDiag, then the map is the identity
        if self.idealCompatible(new_ideal):
            return self


        #We initialize the new list of arcs
        arc_list = []
        #We deep copy each each arc of this ArcDiag's arc list into the new_arc_list
        for single_arc in self.arc_list:
            arc_list.append(SingleArc(single_arc.i,single_arc.j,single_arc.L,single_arc.R,self.n))


        #We iterate through the queue of SingleArc until we find one with a subarc in the input ArcDiag

        for i in range(len(arc_list)):


            for ideal_arc in new_ideal.arc_list:
                #We verify whether the current new ideal SingleArc is subarc of the current SingleArc of this ArcDiag
                if ideal_arc.isSubarc(arc_list[i]):
                    curr_single_arc = arc_list.pop(i)
                    for k in range(ideal_arc.i+1,ideal_arc.j):

                        #We iterate through the internal points of the ideal SingleArc and flip the position of this ArcDiag's SingleArc with regard to the current internal point
                        if k in curr_single_arc.L:
                            curr_single_arc.R.append(k)
                            curr_single_arc.L.remove(k)
                        elif k in curr_single_arc.R:
                            curr_single_arc.L.append(k)
                            curr_single_arc.R.remove(k)
                    #We sort back the L and R lists
                    curr_single_arc.L.sort()
                    curr_single_arc.R.sort()
                    #We add the modified SingleArc back to the queue
                    arc_list.append(curr_single_arc)
                    new_ArcDiag = ArcDiag(arc_list,self.n).untangle()
                    return new_ArcDiag.newBijection(new_ideal)
        """
        for vis_arc in visited_arcs:
            print("~~~~~~~~~~~~~~~~~~~~~~~")
            vis_arc.graph().show()
        print("########################################")
        """
        raise Exception("Not compatible but did not find and arc that had a subarc in the ideal")

    """
    bijectionOneFlipCut

    Bijection in which each arc is flip-cut once, then the diagram is untangled and the process is continued until the ArcDiag is compliant with the input ideal

    Return:
        ArcDiag
    """

    def bijectionOneFlipCut(self,new_ideal):
        #If this ArcDiag is already compatible with the input ArcDiag, then the map is the identity
        if self.idealCompatible(new_ideal):
            return self


        #We initialize the list of arcs 
        arc_list = []

        #We iterate through each arc of this ArcDiag
        for single_arc in self.arc_list:
            #We compare it to the arcs of the ideal ArcDiag until we find an arc that is subarc of the currenct arc
            change = False
            for ideal_arc in new_ideal.arc_list:
                if ideal_arc.isSubarc(single_arc):
                    #We add the flipcut version of the current arc to the list of the new ArcDiag
                    arc_list.append(single_arc.flipCut(ideal_arc.i,ideal_arc.j))
                    change = True
                    break
            #If we looked at every arc of the ideal without changing this arc then it needs no flip cut
            if not change:
                arc_list.append(single_arc)

        #At this step, every arc that needed a flip cut has been flip cut at least once. Now we untangle the ArcDiag and continue
        return ArcDiag(arc_list,self.n).untangle().bijectionOneFlipCut(new_ideal)

    """
    pointToPointFlip

    This method returns the ArcDiag obtained from this ArcDiag by flipping the position of every SingleArc around the points in the input list

    Return:
        ArcDiag
    """

    def pointToPointFlip(self,points_to_flip):
        #This transformation takes this ArcDiag, then for every SingleArc, flips the position of the SingleArc around every points to flip
        new_arc_list = []
        for single_arc in self.arc_list:
            new_L_list = []
            new_R_list = []
            #We look at each interior point in the L set of the current SingleArc. If the point is in the input set, it flips to the R set
            for int_point in single_arc.L:
                if int_point in points_to_flip:
                    new_R_list.append(int_point)
                else:
                    new_L_list.append(int_point)

            #Symmetricaly for the current SingleArc's R set
            for int_point in single_arc.R:
                if int_point in points_to_flip:
                    new_L_list.append(int_point)
                else:
                    new_R_list.append(int_point)

            #Both lists get sorted

            new_L_list.sort()
            new_R_list.sort()
            #We add the newly created SingleArc to the list
            new_arc_list.append(SingleArc(single_arc.i,single_arc.j,new_L_list,new_R_list,self.n))

        return ArcDiag(new_arc_list,self.n)

    """
    obviousBijection

    This method returns the ArcDiag obtained from this ArcDiag where each SingleArc was flipped around each point in the input list followed by an untanglement of the ArcDiag

    Return:
        non crossing ArcDiag

    """

    def obviousBijection(self,points_to_flip):
        flipped = self.pointToPointFlip(points_to_flip)
        flipped.graph(arc_color = 'green',pt_highlight=points_to_flip).show(figsize=1.5)
        return flipped.untangle()



    """
    graph

    This method produces a Graph object that is a drawing of this ArcDiag

    return:
        a sage Graph object representing this ArcDiag
    """

    def graph(self,shift=0,ptsize = 30,pt_highlight = [],arc_color = 'red',anchor=0,direction=0):
        #We start with the points. One point for each value {1,...,n}
        point_list = []
        blue_point_list = []
        n = self.n
        for l in range(1,self.n+1):
            if l in pt_highlight:
                blue_point_list.append([l+shift,0])
            else:
                point_list.append([l+shift,0])

        #We initialize the graph object
        if len(point_list)==0:
            graph = point(blue_point_list,color='blue',size=10+ptsize)
        else:
            graph = point(point_list,color = 'black',size = 10 + ptsize)
            if len(blue_point_list)>0:
                graph+=point(blue_point_list,color='blue',size=10+ptsize)

        #If there is an anchor to display, put it in green
        if anchor >0 and anchor <= self.n:
            graph+= point([anchor+shift,0],color='blue',size=10+ptsize)

        if direction == -1:
            #Display an arrow to the left
            graph+=arrow(((self.n+1)/2+shift+2,-4.5),((self.n+1)/2+shift-2,-4.5),color='purple')

        if direction == 1:
            #Display an arrow to the right
            graph+=arrow(((self.n+1)/2+shift-2,-4.5),((self.n+1)/2+shift+2,-4.5),color='purple')


        #We add each arc to the graph idividually

        for single_arc in self.arc_list:

        #Now we add the arc each semicircle at a time
            semi_circles = single_arc.sepArcs()
            if semi_circles[0][2] == 2:
                #In this case, this SingleArc is a line between two adjacent points
                graph += line([(semi_circles[0][0]+shift,0),(semi_circles[0][1]+shift,0)], color = arc_color,thickness = 3)
            else:
                for semi_circle in semi_circles:
                    curr_half = semi_circle[2]
                    #We use the sage arc graph object for the semicircles
                    graph += arc(((semi_circle[0]+semi_circle[1])/2+shift,0),(semi_circle[1]-semi_circle[0])/2, sector = (0+curr_half*pi,pi+curr_half*pi),thickness=3, color = arc_color)


        #Finally, we draw a box around the drawing
        #print(type(n))
        graph += line([(0+shift,n/2+1),(n+1+shift,n/2+1)])
        graph += line([(0+shift,-(n/2+1)),(n+1+shift,-(n/2+1))])
        graph += line([(0+shift,n/2+1),(0+shift,-(n/2+1))])
        graph += line([(n+1+shift,n/2+1),(n+1+shift,-(n/2+1))])
        graph.axes(False)
        return graph

    """
    meetGraph

    This method produces a Graph object that is a drawing of this meet ArcDiag

    return:
        a sage Graph object representing this ArcDiag
    """
    def meetGraph(self,shift=0,ptsize = 30,pt_highlight = []):
        return self.graph(shift=0,ptsize = 30,pt_highlight = [],arc_color='blue')


    """
    graphPartition

    This method produces a Graph object that is a drawing of this ArcDiag with lines partitionning the interior points in intervals

    return:
        a sage Graph object representing this ArcDiag
    """

    def graphPartition(self,partition,shift=0):
        #We start with the points. One point for each value {1,...,n}
        point_list = []
        n = self.n
        for l in range(1,self.n+1):
            point_list.append([l+shift,0])

        #We initialize the graph object
        graph = point(point_list,color = 'black')

        #We add each arc to the graph idividually
        for single_arc in self.arc_list:

        #Now we add the arc each semicircle at a time
            semi_circles = single_arc.sepArcs()
            if semi_circles[0][2] == 2:
                #In this case, this SingleArc is a line between two adjacent points
                graph += line([(semi_circles[0][0]+shift,0),(semi_circles[0][1]+shift,0)], color = 'red',thickness = 3)
            else:
                for semi_circle in semi_circles:
                    curr_half = semi_circle[2]
                    #We use the sage arc graph object for the semicircles
                    graph += arc(((semi_circle[0]+semi_circle[1])/2+shift,0),(semi_circle[1]-semi_circle[0])/2, sector = (0+curr_half*pi,pi+curr_half*pi),thickness=3, color = 'red')


        for i in range(1,len(partition)):
            graph += line([(shift+partition[i][0]-0.5,-0.25),(shift+partition[i][0]-0.5,0.25)],thickness=2,color = 'green')
        graph += line([(shift+n-0.5,-0.25),(shift+n-0.5,0.25)],thickness=2,color = 'black')
        graph += line([(shift+2-0.5,-0.25),(shift+2-0.5,0.25)],thickness=2,color = 'black')

        #Finally, we draw a box around the drawing
        #print(type(n))
        graph += line([(0+shift,n/2+1),(n+1+shift,n/2+1)])
        graph += line([(0+shift,-(n/2+1)),(n+1+shift,-(n/2+1))])
        graph += line([(0+shift,n/2+1),(0+shift,-(n/2+1))])
        graph += line([(n+1+shift,n/2+1),(n+1+shift,-(n/2+1))])
        graph.axes(False)
        return graph




# ## WOQuotient

# In[7]:


#The "WOQuotient" class is an object that represents a weak order lattice quotient
class WOQuotient:
    """
    A WOQuotient is initialized with:
        a int n: of which symmetric group this WOQuotient is a quotient of
        an ArcDiag ideal: the minimal elements of the forcing order ideal corresponding to this WOQuotient
    """

    def __init__(self,ideal,n):
        self.n = n
        self.ideal = ideal 

    #An WOQuotient of dimension n and k arcs in it's ideal has a string representation of the form:
    #"A quotient on the weak order of S_{n} and ideal with minimal arcs: {ideal.arc_list[0]}, {ideal.arc_list[1]}, ..., {ideal.arc_list[k-1]}"
    def __str__(self):
        return f'A quotient on the weak order of S_{self.n} and ideal with minimal arcs: {" ".join(f"{self.ideal.arc_list[i]}" for i in range(len(self.ideal.
                                                                                                                                                  arc_list)))}'

    def __repr__(self):
        return f'A quotient on the weak order of S_{self.n} and ideal with minimal arcs: {" ".join(f"{self.ideal.arc_list[i]}" for i in range(len(self.ideal.
                                                                                                                                                  arc_list)))}'
    #Two WOQuotient(s) which have the same dimension and same ideal are equal
    def __eq__(self,quo):
        if not isinstance(quo, WOQuotient):
            return NotImplemented

        self.ideal.sort(key = lambda x : x.i)
        quo.ideal.sort(key = lambda x: x.i)

        return self.n == quo.n and self.ideal == quo.ideal




    ########################################
    #### Class Methods #####################
    ########################################

    """
    fromArcList

    This method constructs a WOQuotient object from a list of "arc" tuple of the form "[(L,i,j,R),...]" with "n" providing the dimension

    Return:
        WOQuotient of dimension {n} and ideal with arcs {[(L,i,j,R),...]}
    """

    @classmethod
    def fromArcList(cls,arc_list,n):
        #We simply construct the ideal ArcDiag from the list with the fromArcList class method of the ArcDiag class.
        return cls(ArcDiag.fromArcList(arc_list,n),n)

    """
    fromSingleArcList

    This method constructs a WOQuotient object from a list of SingleArc(s) with "n" providing the dimension

    Return:
        WOQuotient of dimension {n} and ideal with arcs {[(L,i,j,R),...]}
    """

    @classmethod
    def fromSingleArcList(cls,single_arc_list,n):
        #We simply construct the ideal ArcDiag from the list with the fromArcList class method of the ArcDiag class.
        return cls(ArcDiag(single_arc_list,n),n)

    ########################################
    #### Methods ###########################
    ########################################

    """
    listOfIdealArcs

    This method returns the list of all arcs of this WOQuotient's ideal in the tuple representation

    Return:
        A list of tuples (each tuple is an arc of the ideal of this WOQuotient)
    """

    def listOfIdealArcs(self):
        #In this WOQuotient's ideal arc_list, all the arcs are SingleArc objects. Now we want them in the form of tuples
        list_of_arcs = []

        for single_arc in self.ideal.arc_list:
            list_of_arcs.append((single_arc.L,single_arc.i,single_arc.j,single_arc.R))

        list_of_arcs.sort(key = lambda x : x[1])

        return list_of_arcs


    """
    isPermValid

    This method takes as input a permutation (not necessarily of size n) and checks whether this permutation has an
    arc that is part of the ideal of this WOQuotient

    Return:
        True -- if the permutation is valid (no arc in the ideal)
        False -- if the permutation is not valid (at least one arc in the ideal)
    """

    def isPermValid(self,perm):
        #We begin by obtaining the start point of the arcs of the permutation
        descents = Permutation(perm).descents()

        #Now that we have the start points of every arc of the permutation, we can compare these arcs with the arcs in
        #this WOQuotient's ideal

        for perm_arc in descents:
            start_point = (perm[perm_arc],perm[perm_arc-1])
            #print(start_point)
            for ideal_arc in self.ideal.arc_list:
                #An arc A of a permutation may only invalidate this permutation if one of the arcs B in the ideal is a
                #subarc of A. The interval of A thus must contain the interval of B
                if start_point[0] <= ideal_arc.i and start_point[1] >= ideal_arc.j:



                    #An arc A of a permutation may only invalidate this permutation if one of the arcs B in the ideal
                    #is such that B's L set is contained in A's L set and B's R set is contained in A's R set

                    #We need to test for every point in the interval of this ideal arc, so we first assume this ideal arc
                    #is subarc of that permutation arc until we find a point that does not match.
                    is_subarc = True
                    for point in ideal_arc.L:
                        if point in perm[perm_arc+1:]:
                            #In this case, a point is in the L set of this ideal arc but to the right of the descent
                            is_subarc = False
                            break
                    #If we already broke the subarc rule, we don't need to test for the R set
                    if is_subarc:
                        for point in ideal_arc.R:
                            if point in perm[:perm_arc-1]:
                                #In this case, a point is in the R set of this ideal arc but to the left of the descent
                                is_subarc = False
                                break

                    #If this point is reached and that permutation arc never broke the assumption that it has
                    #this ideal arc as a subarc, then we can conclude it is true and the permutation is invalid
                    if is_subarc:
                        return False

        #At this point, we have tested all permutation arcs with all ideal arcs. We can conclude that the permutation
        #is valid

        return True




    """
    minElements2

    This method gives a list of all minimal elements representative of each class of this WOQuotient

    [This second version keeps a list of positions not to try. Not always faster in practice.]

    Return:
        list of permutations (NOT a list of perm objects)
    """

    def minElements2(self):
        #We will obtain the minimal elements by recurrence. For a permutation P of size k that doesn't break the rules 
        #of this quotient, we obtain candidates of size k+1 by adding k+1 at every position P. We test if any of the 
        #arc of this new permutation has an arc in the ideal of this WOQuotient (i.e. any arc in the ArcDiag ideal
        #is a subarc of an arc in the new permutation.

        #We will use two lists: curr_perm_list and next_perm_list. At the kth step, curr_perm_list will contain couples (P,i) where:
        #    P is a permutation that agrees with every contraint of this WOQuotient
        #    i is a list of position in P where adding l > k will break a contraint of WOQuotient

        curr_perm_list = [([1],[])]

        for k in range(2,n+1):
            #We initialize the list of valid permutations of size k
            next_perm_list = []

            #We iterate through all valid k-1 permutations
            for curr_perm in curr_perm_list:
                #print("curr_perm:",curr_perm)
                #We obtain the current permutation's condemned positions
                invalid_pos = curr_perm[1]

                #We initialize a list which will be the valid permutation obtained from the current permutation
                child_list = []

                for pos in range(0,k):
                    if pos in invalid_pos:
                        #In this case we already know this position is condemned
                        continue
                    else:
                        #We create a new candidate
                        #print("curr_perm[0]:",curr_perm[0])
                        #print("pos:",pos)
                        new_perm = list(curr_perm[0])
                        new_perm.insert(pos,k)
                        #print("new_perm:",new_perm)

                        #We test this new candidate
                        if self.isPermValid(new_perm):
                            #If this new permutation passes the test, it is added to the list of children along with the position of k for the offsetting of the invalid positions list
                            child_list.append((new_perm,pos))

                        else:
                            #If this new permutation does not pass the test, then we condemn this position
                            invalid_pos.append(pos)

                #At this point, all positions were tested for this permutation. All we need now is to update the
                #invalid positions for every new permutations to take into account the offsetting by the new element.
                for valid_perm in child_list:
                    #We create a copy of the invalid position list
                    invalid_pos_copy = invalid_pos.copy()
                    for l in range(len(invalid_pos_copy)):
                        if invalid_pos_copy[l] > valid_perm[1]:
                            invalid_pos_copy[l] = invalid_pos_copy[l]+1
                    #Once the invalid list is updated for this permutation, we can send it with its list to the list of next permutations
                    next_perm_list.append((valid_perm[0],invalid_pos_copy))

            #At this point, all permutations went through the process of being added k to all valid places. The list of next permutations becomes the list of current permutations
            curr_perm_list = next_perm_list

        #At this point, we have all valid length n permutation.
        return [curr_perm_list[l][0] for l in range(len(curr_perm_list))]

    """
    minElements

    This method gives a list of all minimal elements representative of each class of this WOQuotient. the perm_object parameter determines if the elements are given as list or as Perm objects

    Return:
        list of permutations/ list of perm objects (if perm_object = True)
    """

    def minElements(self,perm_object = False,print_len=False):
        #We will obtain the minimal elements by recurrence. For a permutation P of size k that doesn't break the rules 
        #of this quotient, we obtain candidates of size k+1 by adding k+1 at every position P. We test if any of the 
        #arc of this new permutation has an arc in the ideal of this WOQuotient (i.e. any arc in the ArcDiag ideal
        #is a subarc of an arc in the new permutation.

        #We will use two lists: curr_perm_list and next_perm_list. At the kth step, curr_perm_list will contain:
        #    P a permutation that agrees with every contraint of this WOQuotient

        curr_perm_list = [[1]]

        for k in range(2,self.n+1):
            #print(k)
            #We initialize the list of valid permutations of size k
            next_perm_list = []

            #We iterate through all valid k-1 permutations
            for curr_perm in curr_perm_list:
                #print("curr_perm:",curr_perm)

                for pos in range(0,k):
                    #We create a new candidate
                    #print("curr_perm[0]:",curr_perm[0])
                    #print("pos:",pos)
                    new_perm = list(curr_perm)
                    new_perm.insert(pos,k)
                    #print("new_perm:",new_perm)

                    #We test this new candidate
                    if self.isPermValid(new_perm):
                        #If this new permutation passes the test, it is added to the list of next permutations
                        next_perm_list.append(new_perm)


                #At this point, all positions were tested for this permutation. 

            #At this point, all permutations went through the process of being added k to all valid places. The list of next permutations becomes the list of current permutations
            curr_perm_list = next_perm_list
        #print(curr_perm_list)
        if print_len:
            print(len(curr_perm_list))

        if perm_object:
            return [Perm(elem) for elem in curr_perm_list]

        #At this point, we have all valid length n permutation.
        return curr_perm_list


    """
    poset

    This method returns the poset object corresponding to this WOQuotient

    Return:
        poset object
    """

    def poset(self,print_len=False):
        #The poset will contain as elements the minimal permutations of this WOQuotient as Perm objects.
        return Poset((self.minElements(True,print_len), lambda x, y: x<=y))

    """
    skeleton

    This method returns the 1-skeleton of this WOQuotient's Hasse diagram

    Return:
        graph object

    """

    def skeleton(self):
        #A 1-skeleton of a directed graph is simply it's undirected counterpart
        return self.poset().hasse_diagram().to_undirected()


    """
    vPolynomial

    This method gives the v-polynomial (Marin, Novelli, Pilaud 2025+) of this WOQuotient.

    Return:
        polynomial in Q{x,y}
    """

    def vPolynomial(self):
        #The v-polynomial lives in the integer polynomial ring which is a subring of the rational polynomial ring
        R.<x,y> = InfinitePolynomialRing(QQ)
        #The v-polynomial registers the start points of the arc diagram of each permutation representative of this WOQuotient
        #We start by obtaining the list of these permutations

        #timeit = time.time()
        list_min = self.minElements()
        #print("minElements time:",time.time()-timeit)


        #The polynomial starts null
        p=0
        #Now for each element of this WOQuotient, we obtain the start points of their arc diagrams
        for elem in list_min:
            #Every element is its own monomial which start constant
            monomial = 1
            #We start by obtaining the descents of the element
            descents = Permutation(elem).descents()
            #For every descent of a permutation, we have an arc going from the second value of the descent to the first
            for des_pos in descents:
                #For every element, its monomial is of the form x_i1*x_i2*...*y_j1*y_j2*... such that in this element's
                #arc diagram representation, one arc starts at x_i1, another starts at x_i2, ..., and one arc ends at y_j1,
                #another ends at y_j2, ...,

                #Since every descent defines an arc (and every arc is defined by a descent), for every descent, this monomial
                #will contain the left value of the descent in the y's and the right value of the descent in the x's
                monomial *= x[elem[des_pos]]*y[elem[des_pos-1]]


            #We went through all descents of this element so its monomial is complete and added to the polynomial
            p+=monomial
            """
            if monomial == x[1]*x[2]*x[3]*y[3]*y[4]*y[5]:
                Perm(elem).graph().show()
            """
        #print(p.coefficient(x[1]*x[2]*y[5]*y[6]))

        #We went through every element of this WOQuotient and thus the polynomial is complete
        return p


    """
    flipClasses

    This method separates the SingleArc(s) of this WOQuotient's ideal ArcDiag into "flip classes" (Marin, Novelli, Pilaud 2025+)

    Return:
        list of (list of SingleArc(s), interval of the union of interior points of the SingleArcs of the class)
    """

    def flipClasses(self):
        #Two SingleArc(s) are in "flip relation" if:
        #    1. they share interior points
        #The partition of the SingleArc(s) into "flip classes" is the transitive closure of this "flip relation"

        #Initialization of the list to be returned
        class_list = []

        #The list of SingleArc(s) is obtained
        list_of_arcs = self.ideal.arc_list

        #We iterate through the SingleArc(s). For each iteration, there are three cases:
        #    1. It is not linked to any previous SingleArc -- A new "flip class" for this SingleArc is created 
        #    2. It is linked to previous SingleArc.(s) of a single "flip class" -- It is added to the "flip class"
        #    3. It is linked to previous SingleArc.(s) of at least 2 "flip classes" -- All of these "flip classes" are merged together and the current SingleArc is added to this "flip class"

        for curr_single_arc in list_of_arcs:
            #We iterate through the classes of class_list to see in which flip class the current SingleArc belongs
            #We keep track of the flip classes this SingleArc is linked to in a list
            curr_class_index_list = []
            for index in range(len(class_list)):
                #We iterate through the SingleArc(s) of this class until we find a SingleArc directly linked to the current SingleArc
                for arc2 in class_list[index]:
                    #We verify whether the two SingleArc(s) are directly linked with the isLinked method from the SingleArc class
                    if curr_single_arc.isLinked(arc2):
                        #In the case where they are directly linked, we add this class' index to the list
                        curr_class_index_list.append(index)
                        break

            #At this point we have gathered all the class indices to which the current SingleArc belongs.
            #In the case where this SingleArc was linked to no previous SingleArc, we create a new "flip class"
            if len(curr_class_index_list) == 0:
                curr_class_index_list.append(curr_single_arc)
                class_list.append(curr_class_index_list)
            #In the case where this SingleArc is linked to a single "flip class" we simply add it to this class
            elif len(curr_class_index_list) == 1:
                class_list[curr_class_index_list[0]].append(curr_single_arc)
            #In the case where the current SingleArc is linked to many "flip classes" we merge them together
            else:
                new_class = []
                while len(curr_class_index_list)>0:
                    new_class+=class_list.pop(curr_class_index_list.pop())
                new_class.append(curr_single_arc)
                class_list.append(new_class)

        #At this point every SingleArc has been checked
        #Now we want to gather for each class, what interval does it span

        classes_with_interval = []
        for arc_class in class_list:
            #We initiate the interval setting the min to be the maximum possible and the max to be the minimum possible
            curr_interval = [self.n,1]
            for single_arc in arc_class:
                #If we find a SingleArc that widens the interval by any direction (or both), we extend the interval [[Badly Written Comment]]
                if single_arc.i+1<curr_interval[0]:
                    curr_interval[0] = single_arc.i+1
                if single_arc.j-1>curr_interval[1]:
                    curr_interval[1] = single_arc.j-1
            #We have gone through every arc and we thus have the correct interval for this class
            classes_with_interval.append([arc_class,[k for k in range(curr_interval[0],curr_interval[1]+1)]])

        return classes_with_interval


    """
    quotientIntClass

    This method separates the interior points {2,...,n-1} of the ideal graph of this WOQuotient into intervals that determines the flip classes of this WOQuotient

    Return:
        Interval partition of {2,...,n-1}
    """

    def quotientIntClass(self):
        #We order the arcs of this WOQuotient's ideal
        self.ideal.sortedArcListStart()

        #We intialize the partition
        partition = []
        #We initialize the first interval
        curr_interval = [2]

        #We also initialize the current position we are checking
        curr_pos = 2
        #for ideal_arc in self.ideal.arc_list:
            #print(ideal_arc)

        #We iterate through all the arcs of the ideal
        for ideal_arc in self.ideal.arc_list:
            #print(ideal_arc)
            #print(curr_pos)
            #First, if the current arc starts at or after the current point, then the current interval is done
            if ideal_arc.i >= curr_pos:
                #print("reached by ", ideal_arc)
                if len(curr_interval) > 0:
                    partition.append(curr_interval)

                #And now each point between the current point and the first interior point of the current ideal_arc are added as singletons
                for i in range(curr_pos+1,ideal_arc.i+1):
                    partition.append([i])

                #Finally, we update the new current interval to be the interior points of the current ideal_arc and update the current position to be the last position of the current interval
                curr_interval = []
                for i in range(ideal_arc.i+1,ideal_arc.j):
                    curr_interval.append(i)
                curr_pos = ideal_arc.j-1

            #Now if the current arc shares interior points with the current interval, we add it's interior points to the interval if necessary
            elif ideal_arc.j-1>curr_pos:
                for i in range(curr_pos+1,ideal_arc.j):
                    curr_interval.append(i)
                #We update the current position
                curr_pos = ideal_arc.j-1
        #Once we have gone through all of the arcs, we add the last interval created
        if len(curr_interval) > 0:
            partition.append(curr_interval)
        #Now we add all the last points which are not interior point of any arc in the ideal
        for i in range(curr_pos+1,self.n):
            #print("reached")
            partition.append([i])
        #print(partition)

        return partition







    """
    quotientFlipClass

    This method lists all lattice quotients of the weak order that are in the same class of JC flip (Marin, Novelli, Pilaud 2025+)
    than this WOQuotient

    Return:
        list of tuples (WOQuotient, list giving the interior points partition
    """

    def quotientFlipClass(self,all_flips=False):
        #A WOQuotient in the same class of JC flip as this WOQuotient is obtained from this WOQuotient by flipping a subset
        #of the SingleArc(s) of the ideal and flipping all SingleArc(s) sharing a "flip class" with at least one SingleArc

        #We begin by grouping this WOQuotient ideal ArcDiag's SingleArc(s) by "flip class"
        #flip_classes is thus a list of list of SingleArcs
        flip_classes = self.flipClasses()
        """
        print("classes:")
        for classes in flip_classes:
            print("********************")
            for single_arc in classes:
                print(single_arc)
        """

        #The list of WOQuotient is initialized
        list_of_quotients = [(self,[])]

        #The other elements are obtained by flipping a subset of the flip classes of this WOQuotient. Some flip classes
        #are symmetric, though, and we have to remove them first

        #We initialize a new list that will contain only the non symmetric flip classes and one that will contain only
        #the symmetric flip_classes
        non_sym_flip_cls = []
        sym_flip_cls = []

        for i in range(len(flip_classes)):
            #We determine whether a certain flip class is symmetric with the isListSymmetric class method of the
            #SingleArc class
            if SingleArc.isListSymmetric(flip_classes[i][0]):
                #Here as these classes are symmetric we don't need to remember their interior points so we throw away the interval
                sym_flip_cls.append(flip_classes[i][0])
            else:
                non_sym_flip_cls.append(flip_classes[i])
        #The symmetric flip classes won't get flipped and their arcs can be grouped together
        sym_single_arcs = []
        for single_arcs in sym_flip_cls:
            sym_single_arcs += single_arcs

        """
        #print("non sym:")
        print("sym:")
        for sym_class in sym_flip_cls:
            print("****************")
            for single_arc in sym_class:
                print(single_arc)
        """

        #print(sym_single_arcs)


        #We now obtain every WOQuotient(s) in this WOQuotient's JC flip class
        if len(non_sym_flip_cls)>0:
            if all_flips:
                for l in range(1,2**(len(non_sym_flip_cls)-1)):
                    #We initialize it's list of arcs with a copy of the arcs that wont get flipped
                    #We also keep track of all the interior points around which a flip occured
                    flipped_int_points = []
                    arc_list = list(sym_single_arcs)
                    #print(arc_list)
                    #We choose a subset to flip
                    bin_l = format(l,f'0{len(non_sym_flip_cls)}b')
                    counter = 0
                    for c in bin_l:
                        if c == '1':
                            #In this case, this flip class is chosen so we flip all of its SingleArc(s)
                            flipped_int_points+=non_sym_flip_cls[counter][1]

                            for single_arc in non_sym_flip_cls[counter][0]:
                                arc_list.append(single_arc.flipped())
                        else:
                            #This flip class was not chosen so we don't flip its SingleArc(s)
                            for single_arc in non_sym_flip_cls[counter][0]:
                                arc_list.append(single_arc)
                        counter+=1
                    #At this point, we have all the SingleArcs of the ideal of a new WOQuotient in the same JC flip class as
                    #this WOQuotient, we can create that new WOQuotient and add it to the list along with a sorted version of the list of interior points flipped to obtain this WOQuotient
                    flipped_int_points.sort()
                    list_of_quotients.append((WOQuotient(ArcDiag(arc_list,self.n),self.n),flipped_int_points))
            else:
                #print("len(non_sym_flip_cls)",len(non_sym_flip_cls))
                for l in range(len(non_sym_flip_cls)):
                    arc_list = list(sym_single_arcs)
                    flipped_int_points = []
                    flipped_int_points+=non_sym_flip_cls[l][1]
                    for j in range(len(non_sym_flip_cls)):
                        if j == l:
                            for single_arc in non_sym_flip_cls[j][0]:
                                arc_list.append(single_arc.flipped())
                        else:
                            for single_arc in non_sym_flip_cls[j][0]:
                                arc_list.append(single_arc)
                    list_of_quotients.append((WOQuotient(ArcDiag(arc_list,self.n),self.n),flipped_int_points))


        return list_of_quotients

    """
    piDown [NEED TO BE FIXED]

    Input: 
        list of permutation

    This method returns a list of Perm objects where the i-th element of the list corresponds to the minimal element of the congruence for this WOQuotient on the weak order of the i-th permutation in the input list

    Return: 
        list of Perm objects
    """

    def piDown(self,input_perm,poset_of_JI,arc_diag_return = False,verbose=False):
        #We obtain the non-crossing arc diagram of the input Perm
        og_arc_diag = input_perm.nonCrossingArcDiag()
        #We obtain the ideal of the poset of MI corresponding to the input Perm
        #print(poset_of_MI.list())
        JI_filter = list(poset_of_JI.order_ideal(og_arc_diag.arc_list))
        #We remove from the current lower ideal the SingleArc(s) contracted by this WOQuotient to obtain X (as defined in (Albertin 2022 p.80 https://theses.hal.science/tel-04212305)
        X = []
        if verbose:
            print(input_perm)
        for filter_arc in JI_filter:
            if verbose:
                print(filter_arc)
            is_in_X = True
            for ideal_arc in self.ideal.arc_list:
                #print("##")
                #print(ideal_arc)
                #print(filter_arc)
                if ideal_arc.isSubarc(filter_arc):
                    is_in_X=False
                    break
            if is_in_X:
                X.append(filter_arc)
        Y = []
        for elt in X:
            #Y is defined as in (Albertin 2022 p.80 https://theses.hal.science/tel-04212305)
            for p in range(elt.i+1,elt.j):
                    #We obtain the sets L \cap ]i,p[, R \cap ]i,p[, L \cap ]p,j[ and R \cap ]p,j[
                    L_p1 = []
                    L_p2 = []
                    R_p1 = []
                    R_p2 = []
                    for val in elt.L:
                        if val < p:
                            L_p1.append(val)
                        elif val > p:
                            L_p2.append(val)
                    for val in elt.R:
                        if val < p:
                            R_p1.append(val)
                        elif val>p:
                            R_p2.append(val)
                        #print("L_p1 = ",L_p1)
                        #print("R_p1 = ",R_p1)
                        #print("L_p2 = ",L_p2)
                        #print("R_p2 = ",R_p2)

                    if SingleArc(elt.i,p,L_p1,R_p1,self.n) in X and SingleArc(p,elt.j,L_p2,R_p2,self.n) in X:
                        #We verify if the current SingleArc is in Y. If so it is removed from the ideal X/Y
                        #print("reached")
                        Y.append(elt)
                        break

        X_minus_Y = [elt for elt in X if elt not in Y]
        #We obtain the minimal elements of the filter of the poset of MI by X\Y
        final_lower_ideal = list(poset_of_JI.order_ideal(X_minus_Y))
        new_filter = list(poset_of_JI.order_ideal_generators(final_lower_ideal))

        #From this new filter, we obtain the ArcDiag
        final_arc_diag = ArcDiag(new_filter,self.n)

        if arc_diag_return:
            return final_arc_diag

        return final_arc_diag.perm()

    """
    piUp [NEED TO BE FIXED]

    Input: 
        Perm perm
        poset poset_of_MI

    This method returns the Perm object corresponding to the pi up of the input Perm object in this WOQuotient

    Return: 
        Perm object if arc_diag_return = False
        ArcDiag if arc_diag_return = True
    """

    def piUp(self,input_perm,poset_of_MI,arc_diag_return = False,verbose = False):
        #We obtain the non-crossing arc diagram of the input Perm
        #print(poset_of_MI.list())
        og_arc_diag = input_perm.meetNonCrossing()
        if verbose:
            print(input_perm)
            print("*(((((((((((((((((((((((((((((((((((((((((((((((((")
            for single_arc in poset_of_MI.list():
                print(single_arc)
            print(")))))))))))))))))))))))))))))))))))))))))))))))))*")
            for single_arc in og_arc_diag.arc_list:
                print(single_arc)
                print(single_arc in poset_of_MI.list())
        #We obtain the ideal of the poset of MI corresponding to the input Perm
        #print("papa")

        #print(poset_of_MI.list())

        MI_filter = list(poset_of_MI.order_ideal(og_arc_diag.arc_list))
        #We remove from the current lower ideal the SingleArc(s) contracted by this WOQuotient to obtain X (as defined in (Albertin 2022 p.80 https://theses.hal.science/tel-04212305)
        X = []
        #print(input_perm)
        for filter_arc in MI_filter:
            #print(filter_arc)
            is_in_X = True
            for ideal_arc in self.ideal.arc_list:
                #print("##")
                #print(ideal_arc)
                #print(filter_arc)
                if ideal_arc.isSubarc(filter_arc):
                    is_in_X=False
                    break
            if is_in_X:
                X.append(filter_arc)
        Y = []
        for elt in X:
            #Y is defined as in (Albertin 2022 p.80 https://theses.hal.science/tel-04212305)
            for p in range(elt.i+1,elt.j):
                    #We obtain the sets L \cap ]i,p[, R \cap ]i,p[, L \cap ]p,j[ and R \cap ]p,j[
                    L_p1 = []
                    L_p2 = []
                    R_p1 = []
                    R_p2 = []
                    for val in elt.L:
                        if val < p:
                            L_p1.append(val)
                        elif val > p:
                            L_p2.append(val)
                    for val in elt.R:
                        if val < p:
                            R_p1.append(val)
                        elif val>p:
                            R_p2.append(val)
                        #print("L_p1 = ",L_p1)
                        #print("R_p1 = ",R_p1)
                        #print("L_p2 = ",L_p2)
                        #print("R_p2 = ",R_p2)

                    if SingleArc(elt.i,p,L_p1,R_p1,self.n) in X and SingleArc(p,elt.j,L_p2,R_p2,self.n) in X:
                        #We verify if the current SingleArc is in Y. If so it is removed from the ideal X/Y
                        #print("reached")
                        Y.append(elt)
                        break

        X_minus_Y = [elt for elt in X if elt not in Y]
        #We obtain the minimal elements of the filter of the poset of MI by X\Y
        final_lower_ideal = list(poset_of_MI.order_ideal(X_minus_Y))
        new_filter = list(poset_of_MI.order_ideal_generators(final_lower_ideal))

        #From this new filter, we obtain the ArcDiag
        final_arc_diag = ArcDiag(new_filter,self.n)

        if arc_diag_return:
            return final_arc_diag

        return final_arc_diag.meetPerm()








    ##
    #Graph Methods
    ##



    """
    idealGraph

    This method produces a Graph object that displays the arc diagram of the minimal elements of the ideal of this
    WOQuotient.

    Return:
        Graph object of an arc diagram of this WOQuotient's ideal.
    """

    def idealGraph(self,shift=0,flip_classes = False):
        if flip_classes:
            return self.ideal.graphPartition(self.quotientIntClass(),shift)
        else:
            return self.ideal.graph(shift)

    def elementGraphSCAB(self,element,JI_poset,MI_poset,shift=0):
        #JI_poset = SingleArc.posetOfJI(self.n)
        #MI_poset = SingleArc.posetOfMI(self.n)

        graph = self.piDown(element,JI_poset,True).graph() + self.piUp(element,MI_poset,True).meetGraph()
        graph.axes(False)
        return graph


