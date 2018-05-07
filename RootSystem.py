#!/usr/bin/env python

def intToCoefStr(i,var, isFirst=False):
    """
    Helper function for correctly rendering integer linear coefficients into strings.
    Used in root system to provide text and TeX representations of roots in terms of a basis.
    E.g.
    """
    if i>1: 
        if isFirst: 
            return str(i)+var
        else:
            return "+ "+str(i)+var
    if i==1:
        if isFirst: 
            return var
        else:
            return "+ "+var
    if i==0: return ""
    if i==-1: return "- "+var
    if i<0: return "- "+str(-i)+var
    return "ERROR"

class Matrix:
    """
    Class representing matrices with integer entries; in particular, for modeling the GCM
    """
    def __init__(self,l):
        if isinstance(l, list):
            for i in range(len(l)):
                n=None
                if isinstance(l[i],list):
                    if n==None:
                        n=len(l[i])
                    elif n!=len(l[i]):
                        raise TypeError("Matrix must be initialized by a list of lists of the same length")
            self.rowDim=len(l)
            self.colDim=len(l[0])
            self.lEntries=l
            self.dEntries={}
            for i in range(self.rowDim):
                for j in range(self.colDim):
                    self.dEntries[(i,j)]=l[i][j]
        else:
            raise TypeError("Matrix must be initialized by a list of lists of the same length")

    @classmethod
    def strToMatrix(self,l):
        pass
    
        
    def getRow(self,i):
        """Returns a vector copy of row i"""
        rowEntries=[]
        #:Blank list to capture row entries
        for j in range(self.colDim):
            rowEntries=rowEntries+[self.dEntries[(i,j)]]
        return Vector(rowEntries[i])

    def getCol(self,j):
        colEntries=[]
        for i in range(self.rowDim):
            colEntries=colEntries+[self.dEntries[(i,j)]]
        return Vector(colEntries)
        
    def __getitem__(self,k):
        if isinstance(k,tuple) and len(k)==2:
            k0,k1=k
            newentries=[]
            if (isinstance(k0,slice)) and (isinstance(k1,slice)):
                for r in self.lEntries[k0]:
                    newrow=[]
                    for e in r[k1]:
                        newrow=newrow+[e]
                    newentries=newentries+[newrow]
            if (isinstance(k0,slice)) and (isinstance(k1,int)):
                for r in self.lEntries[k0]:
                    newentries=newentries+[r[k1]]
            if (isinstance(k0,int)) and (isinstance(k1,slice)):
                r=self.lEntries[k0]
                for e in r[k1]:
                    newentries=newentries+[e]
            if isinstance(k0,int) and isinstance(k1,int):
                return self.dEntries[(k0,k1)]
            if isinstance(newentries[0],list):
                return Matrix(newentries)
            else:
                return Vector(newentries)


    def __mul__(self,other):
        newentries=[]
        if isinstance(other, int):
            for l in self.lEntries:
                newentries=newentries+[[other*x for x in l]]
            return Matrix(newentries) 
        if isinstance(other, Vector):
            if self.colDim!=len(other):
                raise TypeError("Dimension mismatch")
            for i in range(self.rowDim):
                newentries=newentries+[self.getRow(i)*other]
            return Vector(newentries)
        else:
            return NotImplemented

    def __rmul__(self,other):
        newentries=[]
        if isinstance(other, int):
            for l in self.lEntries:
                newentries=newentries+[[other*x for x in l]]
            return Matrix(newentries) 
        if isinstance(other, Vector):
            if self.rowDim!=len(other):
                raise TypeError("Dimension mismatch")
            for i in range(self.colDim):
                newentries=newentries+[self.getCol(i)*other]
        else:
            return NotImplemented
        return Vector(newentries)

    def __str__(self):
        return str(self.lEntries)
    
    def texStr(self):
        st="\\begin{bmatrix}\n"
        for i in self.lEntries:
            for j in i:
                st+=str(j)+" & "
            st= st[:-2]+ "\\\\\n"
        st=st[:-3]+"\n\\end{bmatrix}"
        return st

    def __eq__(self,other):
        if isinstance(other,Matrix): return False
        return self.dEntries==other.dEntries


class Vector:
    
    def __init__(self,l):
        self.entries=l
    @classmethod
    def stdBasVec(self,i,n):
        return Vector(([0]*(i))+[1]+([0]*(n-i-1)))
    def sum(self):
        s=0
        for i in self.entries:
            s+=i
        return s
    def isNonNegative(self):
        n=len(self)
        for i in range(n):
            if (Vector.stdBasVec(i,n)*self)<0:
                return False
        return True
        
    def __getitem__(self,k):
        return self.entries[k]
    
    def __len__(self):
        return len(self.entries)
    
    def __add__(self,other):
        if isinstance(other, Vector):
            if len(self)!=len(other):
                raise TypeError("Can only add vectors of same length")
            newentries=[]
            for i in range(len(self)):
                newentries=newentries+[self[i]+other[i]]
        else:
            return NotImplemented
        return Vector(newentries)
    
    def __sub__(self,other):
        if isinstance(other, Vector):
            if len(self)!=len(other):
                raise TypeError("Can only subtract vectors of same length")
            newentries=[]
            for i in range(len(self)):
                newentries=newentries+[self[i]-other[i]]
        else:
            return NotImplemented
        return Vector(newentries)
    
    def __mul__(self,other):
        if isinstance(other, int):
            newentries=[]
            for i in range(len(self)):
                newentries=newentries+[self[i]*other]
            return Vector(newentries)
        elif isinstance(other, Vector):
            if len(self)!=len(other):
                raise TypeError("Can only multiply vectors of same length")
            prod=0
            for i in range(len(self)):
                prod=prod+self[i]*other[i]
            return prod
        else:
            return NotImplemented
        return len(self.entries)

    def __rmul__(self,other):
        return self*other
    
    def __neg__(self):
        return (-1)*self
    
    def __eq__(self,other):
        return self.entries==other.entries
    
    def __ne__(self,other):
        return not self==other
    
    def __ge__(self,other):
        if self.sum()>other.sum(): return True
        elif self.sum()==other.sum() and self.entries>=other.entries: return True
        else: return False
    
    def __gt__(self,other):
        return self>=other and self!=other
    
    def __le__(self,other):
        return other>=self
    
    def __lt__(self,other):
        return other>self
    
    def __hash__(self):
        return hash(tuple(self.entries))
    
    def __str__(self):
        return str(self.entries)


class RootSystem:
    """
    Models root systems of finite type Lie superalgebras
    """
    def __init__(self, matrix,p):
        """
        Root system is initialized with a symmetrized GCM DA, where A is the GCM and D is the
        diagonal matrix with coprime positive integer diagonal entries.
        """
        self.GCM=matrix
        self.rank=len(self.GCM.lEntries)
        self.roots=[]
        self.rootsByHeight={}
        self.parities={}
        self.sign={}
        self.alpha={}
        self.typestr=None
        for i in range(self.rank):
            self.alpha[i]=Vector.stdBasVec(i,self.rank)
            self.parities[i]=p[i]
            length=self.GCM[i,i]
            if length>=0: self.sign[i]=1
            if length<0: self.sign[i]=-1
            self.roots+=[self.alpha[i]]
        self.generateRoots()

    def reduce(self):
        nRR=self.roots[:]
        for r in nRR:
            d=2*r
            if d in nRR: 
                self.roots.remove(d)
                h=d.sum()
                try:
                    self.rootsByHeight[h].remove(d)
                except KeyError:
                    print(nRR)
                    print(self.rootsByHeight)

    @classmethod
    def buildByType(self, t, p, sign=1, giveReduced=False):
        badInput=False
        typestr=t+str(p)+{1:"+",-1:"-"}[sign]
        if not isinstance(p,list):
            badInput=True
        else:
            if t=="Da":
                pass
            else:
                for i in p:
                    if i!=0 and i!=1: 
                        badInput=True
        if badInput:
            raise TypeError("Second parameter must be a list of 0's and 1's")
        if t=='A':
            rank=len(p)
            m=[[] for i in range(rank)]
            rs=sign
            for i in range(rank):
                if p[i]==0:
                    m[i]+=[rs*2]
                    for j in range(i+1,rank):
                        if j==i+1: 
                            m[i]+=[rs*(-1)]
                            m[j]+=[rs*(-1)]
                        else:
                            m[i]+=[0]
                            m[j]+=[0]
                if p[i]==1:
                    m[i]+=[0]
                    rs=-rs
                    for j in range(i+1,rank):
                        if j==i+1: 
                            m[i]+=[rs*(-1)]
                            m[j]+=[rs*(-1)]
                        else:
                            m[i]+=[0]
                            m[j]+=[0]
            RS=RootSystem(Matrix(m),p)
        elif t=='B':
            rank=len(p)
            m=[[] for i in range(rank)]
            rs=sign
            for i in range(rank-1):
                if p[i]==0:
                    m[i]+=[rs*4]
                    for j in range(i+1,rank):
                        if j==i+1: 
                            m[i]+=[rs*(-2)]
                            m[j]+=[rs*(-2)]
                        else:
                            m[i]+=[0]
                            m[j]+=[0]
                if p[i]==1:
                    m[i]+=[0]
                    rs=-rs
                    for j in range(i+1,rank):
                        if j==i+1: 
                            m[i]+=[rs*(-2)]
                            m[j]+=[rs*(-2)]
                        else:
                            m[i]+=[0]
                            m[j]+=[0]
            m[rank-1]+=[rs*2]
            RS=RootSystem(Matrix(m),p)
        elif t=='C':
            rank=len(p)+1
            m=[[] for i in range(rank)]
            rs=sign
            for i in range(rank-1):
                if p[i]==0:
                    m[i]+=[rs*2]
                if p[i]==1:
                    m[i]+=[0]
                    rs=-rs
                for j in range(i+1,rank):
                    if j==i+1: 
                        if j==rank-1:
                            m[i]+=[rs*(-2)]
                            m[j]+=[rs*(-2)]
                        else:
                            m[i]+=[rs*(-1)]
                            m[j]+=[rs*(-1)]
                    else:
                        m[i]+=[0]
                        m[j]+=[0]
            m[rank-1]+=[rs*4]
            RS=RootSystem(Matrix(m),p+[0])
        elif t=='D':
            rank=len(p)+1
            p=p+[p[-1]]
            m=[[] for i in range(rank)]
            rs=sign
            for i in range(rank-3):
                if p[i]==0:
                    m[i]+=[rs*2]
                if p[i]==1:
                    m[i]+=[0]
                    rs=-rs
                for j in range(i+1,rank):
                    if j==i+1: 
                        m[i]+=[rs*(-1)]
                        m[j]+=[rs*(-1)]
                    else:
                        m[i]+=[0]
                        m[j]+=[0]
            ttlr=[]
            stlr=[]
            lr=[]
            if rank>2:
                if p[-3]==0:
                    ttlr=[2*rs,(-1)*rs,(-1)*rs]
                    stlr=[(-1)*rs]
                    lr=[(-1)*rs]
                else:
                    rs=-rs
                    ttlr=[0,(-1)*rs,(-1)*rs]
                    stlr=[(-1)*rs]
                    lr=[(-1)*rs]

            if p[-1]==0:
                stlr+=[2*rs,0]
                lr+=[0,2*rs]
            else:
                stlr+=[0,2*rs]
                lr+=[2*rs,0]

            if rank>2:
                m[rank-3]+=ttlr
            m[rank-2]+=stlr
            m[rank-1]+=lr

            RS=RootSystem(Matrix(m),p)
        elif t=='Da':
            rank=3
            param=p[0]
            p=p[1]
            if p==[1,0,0]:
                m=[[0,-2,-2*param],
                   [-2,4,0],
                   [-2*param,-0,4*param]]
            if p==[0,1,0]:
                m=[[4,-2,0],
                   [-2,0,2*(param+1)],
                   [0,2*(param+1),-4*(param+1)]]

            RS=RootSystem(Matrix(m),p)
        elif t=='F':
            rank=4
            if p==[1,0,0,0]:
                m=[[0,-1,0,0],
                   [-1,2,-2,0],
                   [0,-2,4,-2],
                   [0,0,-2,4]]
            if p==[0,0,0,0]:
                m=[[2,-1,0,0],
                   [-1,2,-2,0],
                   [0,-2,4,-2],
                   [0,0,-2,4]]
            if p==[0,0,1,0]:
                m=[[4, -2, 0, 0], [-2, 4, -2, 0], [0, -2, 0, 1], [0, 0, 1, 0]]

            RS=RootSystem(Matrix(m),p)
        elif t=='G':
            if p==[1,0,0]:
                m=[[0,-1,0],
                   [-1,2,-3],
                   [0,-3,6]]
            if p==[0,0]:
                m=[[2,-3],
                   [-3,6]]

            RS=RootSystem(Matrix(m),p)

        if giveReduced: RS.reduce()
        RS.typestr=typestr
        return RS

    def rootString(self, r1,r2):
        diff=r1-r2
        if (diff in self.roots) or (-diff in self.roots):
            return self.rootString(diff,r2)+1
        else:
            return 0

    def p(self,vec):
        par=0
        for i in range(len(vec)):
            par+=self.parities[i]*vec[i]
        return par%2

    def isRoot(self,vec):
        return (vec in self.roots)

    def pair(self,v1,v2):
        return (v1*self.GCM*v2)

    def s(self,i,vec=None):
        if vec is None:
            if self.parities[i]==0:
                return self
            newsimples={}
            newparities={}
            for j in range(self.rank):
                newsimples[j]=self.s(i,self.alpha[j])
                newparities[j]=self.p(newsimples[j])
            newentries=[]
            for i in range(self.rank):
                newentries+=[[]]
                for j in range(self.rank):
                    newentries[i]+=[self.pair(newsimples[i],newsimples[j])]
            newmatrix=Matrix(newentries)
            return RootSystem(newmatrix, newparities)
        elif isinstance(vec,Vector):
            if self.GCM[i,i]!=0:
                return vec-(2*(self.alpha[i]*self.GCM*vec)//self.GCM[i,i])*self.alpha[i]
            elif self.isRoot(vec):
                if vec==self.alpha[i]:
                    return (-1)*self.alpha[i]
                elif self.pair(self.alpha[i],vec)==0: return vec
                elif self.isRoot(vec-self.alpha[i]):
                    return vec-self.alpha[i]
                else:
                    return vec+self.alpha[i]



    def generateRoots(self, maxHeight=1000):
        """
        """
        debugcount=0
        #print("Generating roots")
        hGR={}
        largestheight=0
        for h in range(0,maxHeight):
            hGR[h]=[]
        for i in range(self.rank):
            hGR[0]+=[self.alpha[i]]
            largestheight=1
            if self.GCM[i,i]!=0 and self.parities[i]==1:
                hGR[1]+=[2*self.alpha[i]]
                largestheight=2
        self.roots=hGR[0]+hGR[1]
        self.fullRoots=hGR[0]
        for h in range(0,maxHeight):
            #print(h)

            if len(hGR[h])==0: 
                #print("Height "+str(h)+": "+str(hGR[h]))
                break
            for r in hGR[h][:]:
                for i in range(self.rank):
                    v=self.s(i,r)
                    s=v.sum()
                    #print("s_"+str(i)+"("+str(r)+")="+str(v))
                    if (not (v in self.roots)) and h<s-1 and s-1<maxHeight:
                        #print(str(v)+" is not in "+str(self))
                        self.roots+=[v] 
                        hGR[s-1]+=[v]
                        if largestheight<s: largestheight=s
        #print(largestheight)
        self.rootsByHeight={}
        for i in range(0,largestheight): 
            self.rootsByHeight[i+1]=hGR[i]

    def __eq__(self,other):
        if not isinstance(other, RootSystem):
            return False
        if self.rank!= other.rank:
            return False
        if self.parities!=other.parities:
            return False
        if self.GCM!=other.GCM:
            return False
        return True

    def __str__(self):
        s="{ "
        initialRoot=True
        for r in self.roots:
            if not initialRoot: 
                s+=", "
            else: 
                initialRoot=False
            initialSummand=True
            for i in range(len(r)):
                if not (initialSummand and r[i]==0):
                    s+=intToCoefStr(r[i], "a_"+str(i), initialSummand)
                    if initialSummand: initialSummand=False
        s+="}"
        return s
    def texStr(self):
        s="\\{ "
        initialRoot=True
        for r in self.roots:
            if not initialRoot: 
                s+=", "
            else: 
                initialRoot=False
            initialSummand=True
            for i in range(len(r)):
                if not (initialSummand and r[i]==0):
                    s+=intToCoefStr(r[i], "\\alpha_"+str(i), initialSummand)
                    if initialSummand: initialSummand=False
        s+="\\}"
        return s

    def dynkinDiagram(self):
        s="\\begin{tikzpicture}[decoration={markings,mark=at position 0.7 with {\\arrow[black]{>}}}]\n"
        initialRoot=True
        minl=0
        for i in range(self.rank):
            li=abs(self.GCM[(i,i)])
            if li!=0:
                if minl==0:
                    minl=li
                else:
                    minl=min(minl,li)
        if minl==0:
            minl=2
        for i in range(self.rank):
            if self.parities[i]==0:
                s+="\\drawevennode[]("+str(i)+",0)("+str(i)+");\n"
            elif self.GCM[(i,i)]==0:
                s+="\\drawisonode[]("+str(i)+",0)("+str(i)+");\n"
            else:
                s+="\\drawoddnode[]("+str(i)+",0)("+str(i)+");\n"
            if i>0:
                for j in range(0,i):
                    N=abs(self.GCM[(i,j)])
                    li=abs(self.GCM[(i,i)])
                    lin=max(li,2)
                    lj=abs(self.GCM[(j,j)])
                    ljn=max(lj,2)
                    
                    if li!=0 and lj!=0:
                        N= (4*(N**2))//(li*lj)
                    elif (li!=0)^(lj!=0):
                        l=max(li,lj)
                        if l==N:
                            N=2
                    if N==0:
                        s+="%no edges "+str(i)+"--"+str(j)+" because |GCM["+str((i,j))+"]|="+str(abs(self.GCM[(i,j)]))+"\n"
                    elif N==1:
                        if j==i-1:
                            s+="\draw ("+str(j)+".east) -- ("+str(i)+".west);\n"
                        else:
                            s+="\draw ("+str(j)+".south) edge[bend right] ("+str(i)+".south);\n"
                    elif N==2:
                        if j==i-1:
                            if (lj!=0 and li>lj) or (li==0 and 2*N>lj):
                                s+="\draw[line width=.3em, black] ("+str(i)+".west) -- ("+str(j)+".east);\n"
                                s+="\draw[line width=.15em, white,postaction={decorate}] ("+str(i)+".west) -- ("+str(j)+".east);\n"
                            elif (li!=0 and lj>li) or (lj==0 and 2*N>li):
                                s+="\draw[line width=.3em, black] ("+str(j)+".east) -- ("+str(i)+".west);\n"
                                s+="\draw[line width=.15em, white,postaction={decorate}] ("+str(j)+".east) -- ("+str(i)+".west);\n"
                            else:
                                s+="\draw[line width=.3em, black] ("+str(j)+".east) -- ("+str(i)+".west);\n"
                                s+="\draw[line width=.15em, white] ("+str(j)+".east) -- ("+str(i)+".west);\n"
                        else:
                            s+="\draw[] ("+str(j)+".south) edge[bend right] ("+str(i)+".south);\n"
                    else:
                        if j==i-1:
                            s+="\draw ("+str(j)+".east) -- ("+str(i)+".west) node[midway, above] {"+str(N)+"};\n"
                        else:
                            s+="\draw ("+str(j)+".south) edge[bend right] ("+str(i)+".south) node[midway, below] {"+str(N)+"};\n"

                    
        s+="\\end{tikzpicture}\n"
        return s

def buildByType(t, p, sign=1, giveReduced=False):
    return RootSystem.buildByType(t,p,sign,giveReduced)

RS=buildByType("G", [0,0], 1)
print(len(RS.roots))
print(RS)
print(RS.GCM)
l=len(RS.roots)
for i in range(l):
    for j in range(i+1,l):
        ri=RS.roots[i]
        rj=RS.roots[j]
        p=RS.pair(ri,rj)
        if p>0 and ((ri+rj) in RS.roots):
            print(str((str(ri),str(rj),p)))
