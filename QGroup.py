#!/usr/bin/env python
from Wordnomial import Word
from Wordnomial import Wordnomial
from Wordnomial import Monomial
from Wordnomial import Polynomial
from RootSystem import Matrix
from RootSystem import Vector
from RootSystem import RootSystem
import math
import colorsys
import os

def get_N_colors(N=5):

    HSV_tuples = [(x*1.0/N, 0.8, 0.8) for x in range(N)]
    hex_out = []
    s=""
    count=0
    for rgb in HSV_tuples:
        rgb = colorsys.hsv_to_rgb(*rgb)
        s+="\\definecolor{color"+str(count)+"}{rgb}{"+str(rgb)[1:-1]+"}\n"
        count+=1
    return s
def gojs_get_N_colors(N=5):

    HSV_tuples = [(x*1.0/N, 0.6, 0.6) for x in range(N)]
    hex_out = []
    s=[]
    count=0
    for rgb in HSV_tuples:
        rgb = colorsys.hsv_to_rgb(*rgb)
        rgb=[int(100*x) for x in rgb]
        s+=["\"rgb("+str(rgb)[1:-1]+")\""]
        count+=1
    return s

def printN(s):
    printN.N+=1
    print(str(printN.N)+". "+str(s))
printN.N=0

#print(get_N_colors())

def compositions(t=2,s=2):
    """knuth 7.2.1.3"""

    q = [0] * (t+1)
    r = None
    q[0] = s
    while True:
        yield q[:]
        if t==0:
            break
        if q[0] == 0:
            if r==t:
                break
            else:
                q[0] = q[r] - 1
                q[r] = 0
                r = r + 1
        else:
            q[0] = q[0] - 1
            r = 1
        q[r] = q[r] + 1

ONE=Monomial(1,0)
ZERO=Monomial(0,0)
Q=Monomial(0,1)


class QuantumGroup:
    
    """
    A QuantumGroup object facilitates computations in the integral form of the half-quantum enveloping superalgebra
    associated to some root data. 

    Args:
        rS (:mod:`RootSystem`): The first parameter is the root data for the quantum group.

        order (:any:`list`, optional): The second parameter is optional and defaults to None. 
            If no order is given, the default order in rS (i.e. the order induced by the generalized Cartan matrix) is taken. Otherwise,
            it should be a permutation of the numbers :math:`1, 2, ..., r`, where r= :py:attr:`rS.rank` in list form,
            indicating how to reorder the default order.

        useCache (:mod:`bool`, optional): The third parameter is optional and defaults to True. 
            When True, enables caching of shuffle product computations to speed up computation.

    Attributes:
        RS (:class:`RootSystem`): Stores the argument rS. The root system is forced to be reduced in the initialization.

        GCM (:any:`dict`): Stores :py:attr:`rS.GCM` reinterpreted as a :any:`dict`.
            Maps a pair of letters (i,j) to the value of the bilinear pairing :math:`(\alpha_i,\alpha_j)`.

        rank (:any:`int`): Stores local copy of :py:attr:`rS.rank`

        parity (:any:`dict`): Enhanced local copy of :py:attr:`rS.parity`

        alphabet (:any:`list`): List of :py:attr:`rS.rank` characters indexing the simple roots of RS.
            By default, generated automatically in the standard lexicographic order starting with 'a'.

        lexicon (:any:`dict`): Dictionary mapping letters to their lexicographic position.
            By default, the standard lexicographic order is assumed; that is, 'a'->1, 'b'->2, and so on.
            If the optional argument order is used, the function :func:`setLex` is used to reorder the alphabet.

        Lyndon (:any:`dict`): Dictionary mapping roots to their corresponding Lyndon words as a :mod:`Word`.

        LynStr (:any:`dict`): Dictionary mapping roots to their corresponding Lyndon words as a :any:`str`.

        coStdFact (:any:`dict`): Dictionary related to the co-standard factorization of Lyndon words which returns a pair (x,y).
            Let :math:`\ell` be a Lyndon word of weight :math:`r`, and suppose :math:`\ell=uv` is the co-standard factorization into Lyndon
            words :math:`u,v` corresponding to the roots :math:`r_u,r_v`.
            If a Lyndon word :math:`\ell` (:mod:`Word`) is used as a key, then (x,y)=(u,v).
            If a root :math:`r` (:mod:`Word`) is used as a key, then (x,y)=(r_u,r_v).

        RV (:any:`dict`): Dictionary mapping roots to their corresponding root vectors (:mod:`Wordnomial`).

        RVmon (:any:`dict`): Dictionary mapping a root to a pair (x (:mod:`Wordnomial`), y (:mod:`Polynomial`)).
            Here, x is the monomial representation of the root vector (that is, a preimage under :math:`\Psi`)
            and y is a normalization factor.

        RVnorm (:any:`dict`): Dictionary mapping a root to the norm (:mod:`Polynomial`) of the corresponding root vector.

        DRV (:any:`dict`): Dictionary mapping roots to their corresponding dual root vectors (:mod:`Wordnomial`).
            This means :math:`RV[r]=DRV[r]*RVnorm[r]`

        CB (:any:`dict`): Dictionary providing storage for computed canonical basis elements.
            A valid key should be a weight (:mod:`Vector`). If the CB has been computed for that weight,
            it should be returned as a :any:`list` of :mod:`Wordnomial` elements.

        CB (:any:`dict`): Dictionary providing storage for computed dual canonical basis elements.
            A valid key should be a weight (:mod:`Vector`). If the DCB has been computed for that weight,
            it should be returned as a :any:`list` of :mod:`Wordnomial` elements.

        useCache (:any:`bool`): Stores the argument useCache.

        cache (:any:`dict`): Dictionary providing storage for computed shuffle products of words, used if useCache is True.
            Suppose we wish to check if we have computed a shuffle product :math:`w1\diamond w2`. 
            The cache is subdivided into sub-caches based on word length to reduce search time, so we first check for the key k1=(len(w1),len(w2)).
            If the key exists, next we check for the pair of words (represented as strings) k2=(w1.word,w2.word). If both keys exist, then
            :math:`w1\diamond w2= cache[k1][k2]`.

    """

    def __init__(self, rS, order=None, useCache=True, useHDCache=False, debug=False):
        rS.reduce()
        self.debug=debug
        self.useCache=useCache
        self.useHDCache=useCache and useHDCache
        self.nullData() #initialize variables to empty
        self.setRootSystem(rS) # Helper method for correctly setting attributes related to the root system
        if order!= None: self.setLex(order)
        self.genRootVectors(debug=self.debug)

    def nullData(self):
        """
        Sets all of the attributes of the QuantumGroup object to empty.
        """
        self.RS=None
        self.GCM=None
        self.rank=None
        self.parity={}
        self.alphabet=[]
        self.lexicon={}
        self.Lyndon={}
        self.LynStr=[]
        self.coStdFact={}
        self.RV={}
        self.RVmon={}
        self.RVnorm={}
        self.DRV={}
        self.CB={}
        self.DCB={}
        self.cache={}
        self.cacheDiff={}

    def loadCache(self,weight1, weight2, debug=False):
        directory="cache/"+self.RS.typestr+"/"

        if weight1>weight2:
            weight1,weight2=weight2,weight1

        if (weight1,weight2) in self.cache:
            # if debug: 
            #     print("Found "+str((str(weight1),str(weight2)))+" in cache keys, here are the entries")
            #     subcache=self.cache[(weight1,weight2)] 
            #     for x in subcache.keys():
            #         print("\t"+str(x)+": "+str(subcache[x]))
            pass #Don't override the local cache if populated
        elif self.useHDCache:
            filename=directory+str(weight1)+" "+str(weight2)+".txt"
            try:
                with open(filename, "r+") as f:
                    if debug: print("Loading from cache "+filename)

                    temp={}
                    for line in f:
                        line=line[:-1]
                        w1,templine=line.split("*",1)
                        w2,value=templine.split("=")
                        uc=Wordnomial.uncache(value)
                        if debug: print("Loaded \""+line+"\" into "+str(uc))
                        temp[(w1,w2)]=uc
                    self.cache[(weight1,weight2)]=temp
                    self.cacheDiff[(weight1,weight2)]={}

            except IOError as err:
                if debug: print(err)
                self.cache[(weight1,weight2)]={} #If no cached data, initialize local cache and differential cache
                self.cacheDiff[(weight1,weight2)]={}
                pass
        else:
            self.cache[(weight1,weight2)]={} #If no cached data, initialize local cache and differential cache
            self.cacheDiff[(weight1,weight2)]={}

    def appendCache(self, debug=False):
        if not useHDCache:
            return None
        directory="cache/"+self.RS.typestr+"/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        for wep in self.cacheDiff:
            wocache=self.cacheDiff[wep]
            weight1,weight2=wep
            filename=directory+str(weight1)+" "+str(weight2)+".txt"
            try:
                with open(filename, "a+") as f:
                    for wop in wocache:
                        w1,w2=wop
                        f.write(w1+"*"+w2+"="+wocache[wop].cacheStr()+"\n")
            except IOError as err:
                    if debug: print(err)
                    pass 


        
    def printCache(self,MAX=1000):
        """
        Prints the cache. If there are more than MAX(=1000 by default) entries in a bin, or bins, stops printing.
        """
        s="Cache:"
        cache=self.cache
        bincount=0
        for i in sorted(cache.keys()):
            bincount+=1
            s+="\n\n\t Weights "+str([str(x) for x in i])+"\n"
            if bincount>=MAX: 
                break
            subcache=cache[i]
            count=0
            for j in sorted(subcache.keys()):
                if count>0: s+=",\n"
                if count>=MAX: 
                    s+="..."
                    break
                count+=1
                s+="\t\t"+str(j)+" : "+str(subcache[j])
        s+="}"
        print(s)

    def setRootSystem(self,RS):
        """
        Sets the root system, along with other related attributes (rank, alphabet, parity, lexicon, GCM) of the quantum group.

        Args:
            RS (:mod:`RootSystem`): The desired root data for the quantum group.
        """
        if RS is None or RS==self.RS: 
            pass
        else:
            self.nullData() # Clear all old data to prepare for new root system
            RS.reduce() #Reduce root system
            self.RS=RS
            self.rank=RS.rank
            self.alphabet=QuantumGroup.alphByRank(RS.rank)
            parity=RS.parities
            lexicon={}
            for i in range(len(self.alphabet)):
                lexicon[self.alphabet[i]]=i #Order the alphabet following the listr
                parity[self.alphabet[i]]=parity[i]
            self.lexicon=lexicon
            self.parity=parity
            tempA=RS.GCM.dEntries
            A={}
            for key in tempA.keys():
                A[(self.alphabet[key[0]], self.alphabet[key[1]])]=tempA[key]
                A[key]=tempA[key]
            self.GCM=A

    def setLex(self,perm):
        """
        Takes a permutation and sets the lexicographic order accordingly. 

        Args:
            perm (:any:`list` of :any:`int`): a permutation of the list [1,2,3,...,self.rank].

        Note:
            Actually, this is a bit more flexible. Any list of distinct numbers will do, with the new lexicographic
            order determined by comparing these values.
        """
        if isinstance(perm,str):
            np=[]
            for l in self.alphabet:
                minunfilled=0
                try:
                    i=perm.index(l)
                    if minunfilled==i:
                        minunfilled+=1
                    np+=[perm.index(l)]
                except ValueError:
                    np+=[minunfilled]
                    minunfilled+=1
            perm=np
        values=[perm[i] for i in range(len(perm))]
        values.sort()
        for i in range(len(perm)):
            self.lexicon[self.alphabet[i]]=values.index(perm[i])

        self.RV={} 
        self.RVmon={}
        self.DRV={}
        self.CB={}
        self.DCB={}
        self.RVnorm={}
        self.Lyndon={}
        self.coStdFact={}
        self.LynStr=[]
        self.genRootVectors() #Reset all of the data related to the lexicographic order and regenerate it

    def lexKey(self,w):
        """
        Translates a word with respect to a non-standard lexicographic ordering into a corresponding word in the standard lexicographic ordering.
        Essentially, provides an isomorphism of totally ordered sets between the words under the non-standard lex order and the words under
        the standard lex order by applying the corresponding permutation of the alphabet to words letter-by-letter.
        For example, with respect to the order b<a<c, the word abc has the same position as bac in the standard lex order.

        Args:
            w (:any:`str`): a word containing only letters in self.alphabet

        Returns:
            (:any:`str`): a word containing only letters in self.alphabet

        Note:
            This is a quick fix to use Python comparison of strings for arbitrary lexicographic orders, but could easily be broken. Should
            be changed to a less unpredictable method.
        """
        if isinstance(w, str):
            t=""
            for i in w:
                t+=self.alphabet[self.lexicon[i]]
            return t
        else:
            return NotImplemented


    @classmethod
    def alphByRank(cls,rank):
        """
        Returns an alphabet with *rank* letters, starting from 'a'. 

        Args:
            rank (:any:`int`): the desired number of letters in the alphabet (positive)

        Returns:
            (:any:`list` of :any:`chr`): list of letters

        Note:
            For alphabets of more than 26 letters, this method is not desirable. However, computation for such a large alphabet would be unwieldy
            in any case, so it is not clear that it is worth providing for that case.
        """
        return list(map(chr, range(97, 97+rank)))

    def setAlphabet(self,lex):
        """
        Sets a particular alphabet with a particular lexicographic order. 

        Args:
            lex (:any:`dict`): a dictionary whose keys are the letters of the desired alphabet and values are the lexicographic values of the letters
        """
        self.alphabet=lex.keys()
        self.setLex(lex)

    def genLyndonWords(self):
        """
        Uses the root data and alphabet data to generate dominant Lyndon words and sets the attributes :attr:`Lyndon`, :attr:`coStdFact`, and :attr:`LynStr`. 
        This is done using the algorithm described in :doc:`/theory`.

        In short, we must compute a dominant Lyndon word for each root. For simple roots, these are simply letters. We then build inductively
        by height of the root (which is the same as length of the words): we produce all potential costandard factorizations of the corresponding 
        Lyndon word using the known shorter Lyndon words (in particular, they must have the right weight), and the largest product of these 
        must be the Lyndon word.
        """
        RS=self.RS
        lex=self.lexicon
        empty=Word.Empty()
        Lyndon={}
        CSF={}
        alph=self.alphabet
        for i in range(len(alph)): #Initialize the letters as the Length 1 dominant Lyndon Words
            Lyndon[RS.alpha[i]]=Word(Monomial(1,0),[alph[i]])
            CSF[alph[i]]=(alph[i],"")
            CSF[Lyndon[RS.alpha[i]]]=(Lyndon[RS.alpha[i]],Word(Monomial(1,0),[])) #Costandard factorization has trivial right factor in this case
            CSF[RS.alpha[i]]=(RS.alpha[i],0*RS.alpha[i])
        hgr=RS.rootsByHeight
        maxheight=max(hgr.keys())
        for i in range(2,maxheight+1):
            for r in hgr[i]:
                bestcandidate=empty
                bestfactor=(empty,empty)
                bestfactorroots=(0*RS.alpha[0],0*RS.alpha[0])
                for j in range(1,i):
                    for r1 in hgr[j]: #Loops through all roots of less height then r for decompositions r=r1+r2 to find potential costandard factorizations
                        r2=r-r1
                        if r2 in hgr[i-j]:
                            w1=Lyndon[r1]
                            w2=Lyndon[r2]
                            w=w1*w2
                            if Word.compare(w1,w2,lex)==-1 and Word.compare(bestcandidate,w,lex)==-1:  
                                #Determine if w1w2 is a costandard factorization 
                                #and the product is larger than the current best candidate
                                bestcandidate=w
                                bestfactor=(w1,w2)
                                bestfactorroots=(r1,r2)
                Lyndon[r]=bestcandidate
                CSF[bestcandidate]=bestfactor
                CSF[bestcandidate.strWord()]=(bestfactor[0].strWord(),bestfactor[1].strWord())
                CSF[r]=bestfactorroots
        self.Lyndon=Lyndon
        self.coStdFact=CSF
        self.LynStr=[self.Lyndon[root].strWord() for root in self.RS.roots]

    def genRootVectors(self, debug=False):
        """
        Uses the root data and Lyndon data to generate root vectors, their norms, and their duals. 
        Sets the attributes :attr:`RV`, :attr:`DRV`, :attr:`RVmon`, and :attr:`RVnorm`. 
        This is done using the algorithm described in :doc:`/theory`.
        """
        self.genLyndonWords()#Start by ensuring the Lyndon data is determined
        Lyndon=self.Lyndon
        CSF=self.coStdFact
        PBW={} #Dict storing the root vectors in shuffle presentation
        dPBW={} #Dict storing dual root vectors
        PBWmon={} #Dict storing root vectors in monomial presentation
        PBWnorm={} #Dict storing root vector norms under standard bilinear form
        hgr=self.RS.rootsByHeight
        for i in sorted(hgr.keys()):
            for r in hgr[i]:
                if i==1: #Initialize the simplest case: height one root vectors are simply the generators
                    PBW[r]=Lyndon[r]
                    PBWmon[r]=(Lyndon[r],1)
                    PBWnorm[r]=1
                else: #Otherwise, the root vector is inductively built out of the root vectors in the costandard factorization
                    r2=CSF[r][1]
                    r1=CSF[r][0]
                    pair1=self.RS.pair(r1,r1)
                    pair2=self.RS.pair(r2,r2)
                    pair3=self.RS.pair(r2,r1)

                    # Determine relative ``lengths'', even if one or both roots are isotropic.
                    #if pair1==0 and pair2!=0:
                    #    l1=abs((2*pair3)//pair2)
                    #    l2=abs(pair2)
                    #elif pair2==0 and pair1!=0:
                    #    l2=abs((2*pair3)//pair1)
                    #    l1=abs(pair1)
                    #else:
                    #    l1,l2=1,1

                    #if pair1>pair2:
                    #    cr=r2
                    #    length=l2
                    #    pair=pair2
                    #else:
                    #    cr=r1
                    #    length=l1
                    #    pair=pair1
                    pair=min(pair1,pair2)

                    for j in range(self.RS.rank):
                        alph=self.RS.alpha[j]
                        ap=self.RS.pair(alph,alph)
                        if abs(ap)==abs(pair):
                            length=max(abs(ap),1)
                            break

                    if self.RS.p(alph)==1 and ap!=0:
                        sign=-1
                    else:
                        sign=1

                    if length%2==0: length=length//2
                    
                    string=self.RS.rootString(r2,r1)
                    if debug:
                        print("\n"+
                            "----------------------------\n"+
                            "Computing PBW vector for root "+str(r) +"\n"+
                            "Lyndon is "+str(Lyndon[r]) +" with CSF "+str((str(Lyndon[r1]),str(Lyndon[r2])))+"\n"+
                            "Additional normalization constant is [" + str(string+1)+"]_{q^"+str(length)+"} = "+str(Polynomial.qInt(string+1,sign,length=length)))

                    PBW[r]=self.bracket((PBW[r2]),(PBW[r1]), shuffle=True).normalize(Polynomial.qInt(string+1,sign,length=length))

                    if debug:
                        print("\n"+
                            "Root Vector is [ RV["+str(r2)+"], RV["+str(r1)+"] ]_q / ["+str(string)+"+1]_{q^"+str(length)+"} = "+str(PBW[r])
                            )

                    PBWmon[r]=(self.bracket((PBWmon[r2][0]),(PBWmon[r1][0]), shuffle=False),PBWmon[r1][1]*PBWmon[r2][1]*Polynomial.qInt(string+1,sign,length=length))
                    #Note for monomial presentation, the bracket product uses concatenation rather than the shuffle product
                    PBWnorm[r]=Wordnomial([])
                    for w in PBWmon[r][0].listRep():
                        coef=(w.coefficient())
                        pair=PBW[r][w.strWord()[::-1]].coefficient()
                        PBWnorm[r]+=coef*pair
                    try:
                        PBWnorm[r]=PBWnorm[r].normalize(PBWmon[r][1])
                    except ArithmeticError:
                        print("PBWnorm computation error:\n PBWmon["+str(r)+"][0]="+str(PBWmon[r][0])
                              +"\n PBWmon["+str(r1)+"][1]="+str(PBWmon[r1][1])
                              +"\n PBWmon["+str(r2)+"][1]="+str(PBWmon[r2][1])
                              +"\n PBWmon["+str(r)+"][1]="+str(PBWmon[r][1])
                              +"\n PBWnorm["+str(r)+"]="+str(PBWnorm[r])
                              +"\n Normalization factor="+str(PBWmon[r][1])
                              +"\n String length="+str(string)
                              +"\n Qint factor="+str(Polynomial.qInt(string+1,sign,length=length))
                              +"\n Lyndon Words="+str([l for l in self.LynStr])
                              )
                        self.printCache()
                        quit()
                    try:
                        PBWnorm[r]=PBWnorm[r].listRep()[0]
                    except IndexError:
                        print("Error: Norm of PBW["+str(r)+"]="+str(PBWmon[r][0])+"=["+str(PBWmon[r2][0])+","+str(PBWmon[r1][0])+"] is zero!")
                        quit()
                    PBWnorm[r]=PBWnorm[r].coefficient()

        for r in PBWnorm:
            if isinstance(PBW[r],Wordnomial):
                try:
                    dPBW[r]=PBW[r].normalize(PBWnorm[r])
                except ArithmeticError:
                        print("dPBW computation error:\n PBWmon["+str(r)+"][0]="+str(PBWmon[r][0])
                              +"\n PBW["+str(r)+"]="+str(PBW[r])
                              +"\n PBWmon["+str(r1)+"][1]="+str(PBWmon[r1][1])
                              +"\n PBW["+str(r1)+"]="+str(PBW[r1])
                              +"\n PBWmon["+str(r2)+"][1]="+str(PBWmon[r2][1])
                              +"\n PBW["+str(r2)+"]="+str(PBW[r2])
                              +"\n PBWmon["+str(r)+"][1]="+str(PBWmon[r][1])
                              +"\n PBWnorm["+str(r)+"]="+str(PBWnorm[r])
                              +"\n Normalization factor="+str(PBWmon[r][1])
                              +"\n String length="+str(string)
                              +"\n Qint factor="+str(Polynomial.qInt(string+1,sign,length=length))
                              +"\n Lyndon Words="+str([l for l in self.LynStr])
                              )
            else: dPBW[r]=PBW[r]
            if debug:
                        print("\n"+
                            "----------------------------\n"+
                            "Computing dual PBW vector for root "+str(r) +"\n"+
                            "Norm is "+str(PBWnorm[r]) +"\n"+
                            "Dual Root Vector is  "+str(dPBW[r])
                            )
           # if debug:
           #     print([str(x) for x in CSF])
           #     print("DEBUG:"
           #       +"\n Lyndon["+str(r)+"]="+str(Lyndon[r])
           #       +"\n PBW["+str(r)+"]="+str(PBW[r])
           #       +"\n CSF["+str(r)+"]="+str(CSF[r][0])
           #       +"\n PBWmon["+str(r)+"][0]="+str(PBWmon[r][0])
           #       +"\n PBWmon["+str(r)+"][0]="+str(PBWmon[r][0])
           #       +"\n PBWmon["+str(r)+"][1]="+str(PBWmon[r][1])
           #       +"\n PBWnorm["+str(r)+"]="+str(PBWnorm[r])
           #       +"\n PBW^*["+str(r)+"]="+str(dPBW[r]))

        self.RV=PBW
        self.DRV=dPBW
        self.RVmon=PBWmon
        self.RVnorm=PBWnorm



    def rootLattice(self,v):
        """
        Takes a weight vector and converts it to an element of the root lattice

        Args:
            w (:any:`str` or :any:`list` or :mod:`Word`): The word whose weight we wish to compute
        """
        s=""
        ent=v.entries
        for x in range(len(v)):
            e=ent[x]
            if e==0: pass
            else:
                if x>0 and len(s)>0: s+="+"
                if e<0:
                    s+="-"
                    e=-e
                if e>1:
                    s+=str(e)
                s+=self.alphabet[x]
        return s



    def weight(self,w):
        """
        Takes a word (either as a string, list, or Word) and computes the corresponding weight vector.

        Args:
            w (:any:`str` or :any:`list` or :mod:`Word`): The word whose weight we wish to compute
        """
        weightvec=[0 for i in range(self.rank)]
        if isinstance(w, str) or isinstance(w, list): word=w
        if isinstance(w, Word):word=w.getWord()
        for l in word:
            weightvec[self.alphabet.index(l)]+=1
        return Vector(weightvec)

    def weightDecomp(self, w):
        """
        Takes a wordnomial, and splits it up into weightwise components

        Args:
            w (:mod:`Wordnomial`): The wordnomial to be decomposed

        Returns:
            (:mod:`dict` of :mod:`Vector`->:mod:`Wordnomial`): dict mapping a weight to the weight component
        """
        #
        l=w.listRep()
        wd={}
        for x in l:
            weight=self.weight(x)
            if weight in wd.keys():
                wd[weight]+=x
            else:
                wd[weight]=x
        return wd

    def weightpair(self,w1,w2):
        """
        Takes two words w1, w2 (either as a string, list, or Word) and compute the value of the bilinear pairing (weight(w1),weight(w2))

        Args:
            w1 (:any:`str` or :any:`list` or :mod:`Word`): The first word in the pair
            w2 (:any:`str` or :any:`list` or :mod:`Word`): The second word in the pair

        Returns:
            (:any:`int`): the value of (weight(w1),weight(w2)) with respect to the root system :attr:`RS`
        """
        return self.RS.pair(self.weight(w1),self.weight(w2))

    def sign(self,w1,w2):
        """
        Takes two words w1, w2 (either as a string, list, or Word) and compute the parity sign (-1)**(p(w1)p(w2))

        Args:
            w1 (:any:`str` or :any:`list` or :mod:`Word`): The first word in the pair
            w2 (:any:`str` or :any:`list` or :mod:`Word`): The second word in the pair

        Returns:
            (:any:`int`): the value of (-1)**(p(w1)p(w2)) with respect to the root system :attr:`RS`
        """
        return (-1)**(self.RS.p(self.weight(w1))*self.RS.p(self.weight(w2)))

    def bar(self,w, debug=False):
        """
        Takes a :mod:`Wordnomial` (or automatically converts to one from an :any:`int`, :mod:`Monomial`, :mod:`Polynomial`, or :mod:`Word`) 
        and computes its image under the bar involution; see :doc:`/theory` for more information.

        Args:
            w (:any:`int`, :mod:`Monomial`, :mod:`Polynomial`, :mod:`Word`, or :mod:`Wordnomial`): The wordnomial
            debug (:any:`bool`, optional, defaults to False): Option for printing debug messages

        Returns:
            (:mod:`Wordnomial`): the image of w under the bar involution
        """
        if debug: print("Debugging bar:")
        if isinstance(w, int):
            if debug: print("Converting integer"+str(w)+" to Word ")
            w=Word(Monomial(w,0),[])
        if isinstance(w, Monomial) or isinstance(w, Polynomial):
            if debug: print("Converting Polynomial"+str(w1)+" to Word ")
            w1=Word(w,[])

        if isinstance(w, Word):
            coef=w.coefficient()
            bcoef=coef.bar()
            length=len(w)
            lr=w.getWord()
            spow=0
            qpow=0
            for i in range(length):
                li=lr[i]
                for j in range(i+1,length):
                    lj=lr[j]
                    qpow-=self.GCM[(li,lj)]
                    spow+=self.parity[li]*self.parity[lj]
            bcoef=Monomial((-1)**spow,qpow)*bcoef
            return Word(bcoef,lr[::-1])

            
        if isinstance(w,Wordnomial):
            new = Wordnomial([])
            for each in w.listRep():
                    new=new+self.bar(each)
            return new
        else:
            return NotImplemented

    def tau(self,w, debug=False):
        """
        Takes a :mod:`Wordnomial` (or automatically converts to one from an :any:`int`, :mod:`Monomial`, :mod:`Polynomial`, or :mod:`Word`) 
        and computes its image under the standard anti-involution; see :doc:`/theory` for more information.

        Args:
            w (:any:`int`, :mod:`Monomial`, :mod:`Polynomial`, :mod:`Word`, or :mod:`Wordnomial`): The wordnomial
            debug (:any:`bool`, optional, defaults to False): Option for printing debug messages

        Returns:
            (:mod:`Wordnomial`): the image of w under the anti-involution
        """
        if debug: print("Debugging bar:")
        if isinstance(w, int):
            if debug: print("Converting integer"+str(w)+" to Word ")
            w=Word(Monomial(w,0),[])
        if isinstance(w, Monomial) or isinstance(w, Polynomial):
            if debug: print("Converting Polynomial"+str(w)+" to Word ")
            w=Word(w,[])

        if isinstance(w, Word):
            return w.revWord()

            
        if isinstance(w,Wordnomial):
            new = Wordnomial([])
            for each in w.listRep():
                    new=new+self.tau(each)
            return new
        else:
            return NotImplemented
        


    def shuffle(self, w1,w2, debug=False):
        """
        Takes two w1,w2 of type :mod:`Word` or :mod:`Wordnomial` (or automatically converts to one from an :any:`int`, :mod:`Monomial`, or :mod:`Polynomial`) 
        and computes the shuffle product of w1 with w2

        Args:
            w1 (:any:`int`, :mod:`Monomial`, :mod:`Polynomial`, :mod:`Word`, or :mod:`Wordnomial`): The first wordnomial factor
            w2 (:any:`int`, :mod:`Monomial`, :mod:`Polynomial`, :mod:`Word`, or :mod:`Wordnomial`): The second wordnomial factor
            debug (:any:`bool`, optional, defaults to False): Option for printing debug messages

        Returns:
            (:mod:`Wordnomial`): the shuffle product of w1 with w2
        """
        if debug: print("Debugging shuffle:")
        if isinstance(w1, int):
            if debug: print("Converting integer"+str(w1)+" to Word ")
            w1=Word(Monomial(w1,0),[])
        if isinstance(w1, Monomial) or isinstance(w1, Polynomial):
            if debug: print("Converting Polynomial"+str(w1)+" to Word ")
            w1=Word(w1,[])
        if isinstance(w2, int):
            if debug: print("Converting integer"+str(w2)+" to Word ")
            w1=Word(Monomial(w2,0),[])
        if isinstance(w2, Monomial) or isinstance(w2, Polynomial):
            if debug: print("Converting Polynomial"+str(w2)+" to Word ")
            w1=Word(w2,[])

        if isinstance(w1, Word) and isinstance(w2, Word):
            if debug: print("Shuffling "+str(w1)+" with "+str(w2))
            if len(w1)==0 or len(w2)==0: return w1*w2
            else:
                c=w1.coefficient()*w2.coefficient()
                w=w1.word[:]
                ow=w2.word[:]
                sw="".join(w)
                sow="".join(ow)
                lenw=len(w)
                lenow=len(ow)
                we=self.weight(sw)
                owe=self.weight(sow)
                if self.useCache:
                    if debug: 
                        print("Searching Cache")
                    rev=we>owe
                    if debug: 
                        print(str(we)+">" +str(owe)+" is " +str(rev))
                        if rev: print("Reversing mult order for retrieving cache")
                    self.loadCache(we,owe)
                    if rev: subcache=self.cache[(owe,we)]
                    else: subcache=self.cache[(we,owe)]
                    keys=subcache.keys()
                    if debug: 
                        print("Keys found in cache: "+str([str(x) for x in keys]))
                    if rev and ((sow[::-1],sw[::-1]) in keys):
                        shf= self.tau(c*subcache[(sow[::-1],sw[::-1])])
                        if debug: print("after tau, retrieved from cache: "+sw+"*"+sow+"="+str(shf))
                        return shf
                    elif (not rev) and (sw,sow) in keys:
                        shf= c*subcache[(sw,sow)]
                        if debug: print("retrieved from cache: "+sw+"*"+sow+"="+str(shf))
                        return shf
                    else:
                        if debug: 
                            print("Not found in cache, computing directly")
                        pass
                l=w[-1]
                ol=ow[-1]
                fshuffle=self.shuffle(Word(Monomial(1,0),w[0:len(w)-1]),Word(Monomial(1,0),ow))
                if debug: 
                    print("fshuffle="+str(fshuffle))
                sshuffle=self.shuffle(Word(Monomial(1,0),w),Word(Monomial(1,0),ow[0:len(ow)-1]))
                if debug: 
                    print("sshuffle="+str(fshuffle))
                newnomial=(fshuffle*Word(Monomial(1,0),[l])+(Monomial(self.sign(w,ol),-self.weightpair(w,ol))*sshuffle)*Word(Monomial(1,0),[ol]))
                if self.useCache: 
                    if rev:
                        temp=self.tau(newnomial)
                        self.cache[(owe,we)][(sow[::-1],sw[::-1])]=temp
                        if self.useHDCache: self.cacheDiff[(owe,we)][(sow[::-1],sw[::-1])]=temp
                    else:
                        self.cache[(we,owe)][(sw,sow)]=newnomial
                        if self.useHDCache: self.cacheDiff[(we,owe)][(sw,sow)]=newnomial
                newnomial=c*newnomial
                if debug: 
                    print("Shuffle resulted in "+str(newnomial))
                return newnomial
            
        if (isinstance(w1,Word) or isinstance(w1,Wordnomial)) and (isinstance(w2,Word) or isinstance(w2,Wordnomial)):
            if isinstance(w1,Word): w1=Wordnomial([w1])
            if isinstance(w2,Word): w2=Wordnomial([w2])
            new = Wordnomial([])
            for each in w1.listRep():
                for other in w2.listRep():
                    new=new+self.shuffle(each,other, debug)
            return new
        else:
            return NotImplemented
    
    def bracket(self, w1, w2,sign=-1, shuffle=True): 


        if isinstance(w1,Word)  and isinstance(w2,Word) :
            if shuffle:
                x=self.shuffle(w1,w2)
                y=self.shuffle(w2,w1)
            else:
                x=w1*w2
                y=w2*w1
            #print("["+str(w1)+","+str(w2)+"] = ("+str(x)+") - "+str((Monomial(self.sign(w1,w2),sign*self.weightpair(w1,w2))))+" (" +str(y)+") = "+str(x+(Monomial(-self.sign(w1,w2),sign*self.weightpair(w1,w2)))*y))
            return x+(Monomial(-self.sign(w1,w2),sign*self.weightpair(w1,w2)))*y
        if (isinstance(w1,Word) or isinstance(w1,Wordnomial)) and (isinstance(w2,Word) or isinstance(w2,Wordnomial)):
            if isinstance(w1,Word): w1=Wordnomial([w1])
            if isinstance(w2,Word): w2=Wordnomial([w2])
            newnomial=Wordnomial([Word(0,[])])
            for each in w1.listRep():
                for other in w2.listRep():
                    newnomial=newnomial+self.bracket(each,other,sign, shuffle)
            return newnomial

    def divpow(self,x,n, sign=1,length=1, debug=False): 
        if n<0: return Wordnomial([Word(Monomial(0,0),[])])
        if n==0: return Wordnomial([Word(Monomial(1,0),[])])
        rec=self.shuffle(x,self.divpow(x,n-1,sign,length))
        if debug: print("divpow debug: rec="+str(rec))
        temp=rec.normalize(Polynomial.qInt(n,sign,length))
        if debug: print("divpow debug: normalizer="+str(Polynomial.qInt(n,sign,length,debug)))
        if debug: print("divpow debug: temp="+str(temp))
        return temp

    def nondivpow(self,x,n): 
        if n<0: return Wordnomial([Word(Monomial(0,0),[])])
        if n==0: return Wordnomial([Word(Monomial(1,0),[])])
        temp=self.shuffle(x,self.nondivpow(x,n-1))  
        return temp

    def PBWpowerdict(self):
        d=dict() 
        words=self.LynStr
        for x in words: d[x]=0
        return d

    def PBWpowerWeight(self,d):
        weight=0*self.RS.alpha[0]
        for k in d.keys():
            weight+=d[k]*self.weight(k)
        return weight

    def PBWpowers(self,weight, index=0, d={}, debug=False): #takes a weight vector, returns a list of dictionaries lyndon->power representing a power vector of a PBW basis elt of that weight
        if index>=len(self.RS.roots):
            return [d]
        root=self.RS.roots[index]
        lr=self.Lyndon[root].strWord()
        d[lr]=0
        valid=weight.isNonNegative()
        l=[]
        while valid:
            #print(str(weight)+"\t"+str(index)+"\t"+str(root)+"\t"+str(d[lr]))
            #print(d)
            nl=self.PBWpowers(weight=(weight-d[lr]*root),index=index+1,d=d.copy(), debug=debug)
            if index==0:
                for nd in nl:
                    iso=False
                    if debug: print(d)
                    for lyn in nd.keys():
                        iso = iso or (self.weightpair(lyn,lyn)==0 and nd[lyn]>1)
                    nw=self.PBWpowerWeight(nd)
                    if debug:
                        print(str(nw)+"=="+str(weight)+"\t"+str(nw==weight)+"\t Iso? " + str(iso))
                    if nw==weight and not iso:
                        l+=[nd]
            else: l+=nl
            d[lr]+=1
            valid=(weight-d[lr]*root).isNonNegative()
            #print(str(weight-d[lr]*root)+"\t"+str(index)+"\t"+str(root)+"\t"+str(d[lr]))
        return l


    def PBWstring(self,d): #takes a dict of pbw powers, produces corresponding dominant word
        word=""
        words=self.LynStr
        words=sorted(words,key=self.lexKey, reverse=True)
        for l in words:
            word=word+l*d[l]
        return word


    def PBWwords(self,weight): #returns a list of dominant words of a fixed weight
        x=self.PBWpowers(weight)
        words=[]
        for d in x:
            words=words+[self.PBWstring(d)]
        return words
        
    def PBWwordpowDict(self, weight): #returns a dictionary of dominant words -> pbw power dicts of a fixed weight
        x=self.PBWpowers(weight)
        wordpow={}
        for d in x:
            s=self.PBWstring(d)
            wordpow[s]=d
        return wordpow

    def PBWel(self, d, dual=False): #takes a dict of pbw powers or a string representing a dominant word, produces corresponding PBW
        if isinstance(d,str):
            wt=self.RS.alpha[0]*0
            for i in d:
                wt+=self.weight(i)
            ddict=self.PBWwordpowDict(wt)
            if d in ddict:
                d=ddict[d]
            else:
                print(d+" is not a dominant word! Returning None.")
                return None
            
        words=self.LynStr
        words=sorted(words,key=self.lexKey)
        pbw=Wordnomial([Word(Monomial(1,0),[])])
        for l in words:
            weight=self.weight(l)
            pair=self.RS.pair(weight,weight)
            sign=(-1)**(self.RS.p(weight))
            length=max(abs(pair)//2,1)
            power=d[l]
            if power>0:
                if pair<0: 
                    if dual:
                        #print((l,str(power),str((((power-1)*(power-2))//2))))
                        pbw=Monomial(sign**(((power-1)*(power-2))//2),(pair*(power*power-power))//4)*self.shuffle(pbw,self.divpow(self.DRV[weight],power, sign=sign, length=length))
                    else:
                        pbw=self.shuffle(pbw,self.nondivpow(self.RV[weight],power))
                else: 
                    if dual:
                        pbw=Monomial(sign**(((power-1)*(power-2))//2),(pair*(power*power-power))//4)*self.shuffle(pbw,self.divpow(self.DRV[weight],power, sign=sign, length=0))
                    else:
                        pbw=self.shuffle(pbw,(self.divpow(self.RV[weight],power, sign=sign, length=length)))
        return pbw

    def MonEl(self, w,div=False): #takes a string, produces corresponding monomial
        mon=Wordnomial([Word(Monomial(1,0),[])])
        lc=0
        ll=w[0]
        for l in w+" ":
            if ll==l:
                lc+=1
            else:
                weight=self.weight(ll)
                length=self.RS.pair(weight,weight)//2
                sign=(-1)**(self.RS.p(weight))
                if div: 
                    mon=self.shuffle(mon,(self.divpow(self.PBWel(ll),lc, sign=sign, length=length)))
                else: 
                    mon=self.shuffle(mon,(self.nondivpow(self.PBWel(ll),lc)))
                ll=l
                lc=1
        return mon




    def PBWexp(self,x, weight=None, string=False, dual=False, debug=False): #x is a wordnomial of weight a alpha_1 + b alpha_2 + c alpha_3
        if debug: print("PBWexp in debug mode")
        if weight is None:
            wd=self.weightDecomp(x)
            weights=wd.keys()
            s=""
            expansion=Wordnomial([])
            initial=True
            for weight in weights:
                if not initial: 
                    s+= " + "
                else:
                    initial=False
                if string: s+=self.PBWexp(wd[weight],weight, string, dual, debug)
                else: 
                    try:
                        pbwe=self.PBWexp(wd[weight],weight, string, dual, debug)
                        expansion+=pbwe
                    except TypeError:
                        raise TypeError("Expansion"+str(pbwe))
            if string: return s
            else: return expansion

        PBWd=self.PBWpowers(weight)
        PBWw=[]
        PBWmult=dict()
        y=x
        for d in PBWd:
            s=self.PBWstring(d)
            PBWw=PBWw+[s]
        PBWw=sorted(PBWw,key=self.lexKey, reverse=True)
        if debug: print(PBWw)
        m=""
        expansion=Wordnomial([])
        init=True
        for w in PBWw:
            if isinstance(y,Word): y=Wordnomial([y])
            wcoef=(y[w]).coefficient()
            pbwel=self.PBWel(w, dual)
            pbwcoef=pbwel[w].coefficient()
            if pbwcoef.isZero(): 
                raise ArithmeticError("Coefficient of "+w+" in PBW("+w+") is zero!?")
            PBWmult[w]=(wcoef/pbwcoef)[0]
            y=y-(PBWmult[w]*pbwel)
        PBWw=sorted(PBWw,key=self.lexKey)
        for w in PBWw:
            pre="PBW"
            if dual: pre+="^*"
            expansion=expansion+Word(PBWmult[w],[pre+":",w])
            if not PBWmult[w].isZero():
                if init: init=not init
                else: m=m+"+ "
                if (PBWmult[w]-Monomial(1,0)).isZero(): m=m+pre+"("+w+")"
                else: m=m+"( "+str(PBWmult[w])+" ) "+pre+"("+w+")"
        if y.isZero(): 
            if string: return m
            else: return expansion
        else: 
            print("PBW expansion Error: "+str(x)+"=("+m+")+"+str(y))
            print(self.PBWel("abb"))
            print(self.RV[Vector([1,2])])
            print(self.MonEl("ba", div=True)-Monomial.Q(2)*self.MonEl("ab", div=True))
            print(self.PBWel("ab"))
            print(self.MonEl("bba", div=True)-Monomial.Q(1)*self.MonEl("bab", div=True)+Monomial.Q(2)*self.MonEl("abb", div=True))
            quit()

    def qDif(self,i,w, left=True):
        if left: ind=-1
        else: ind=0   
        if isinstance(i,int):
            i=self.alphabet[i]
        nw=Wordnomial([])
        if not isinstance(w,Wordnomial): raise ArithmeticError("Not wordnomial")
        for word in w.listRep():
            if not isinstance(word,Word): raise ArithmeticError("Not word: "+str(w))
            s=word.getWord()
            if s!= [] and s[ind]==i:
                if left: nw=nw+Word(word.coefficient(),s[:-1])
                else: nw=nw+Word(word.coefficient(),s[1:])
        return nw

    def bil(self,x,y):
        if isinstance(x,Word) and isinstance(y,Word):
            if x.isEmpty():
                if y.isEmpty():
                    return x.coefficient()*y.coefficient()
                else:
                    return Wordnomial([])
            else:
                xsl=x.strWord()[0]
                return self.bil(self.qDif(xsl,x,left=False),self.qDif(xsl,y,left=True))





    def maxKash(self,i,w):
        if isinstance(i,int):
            i=self.alphabet[i]
        y=w
        ey=self.qDif(i,y)
        n=0
        while not (ey.isZero()):
            n=n+1
            y=ey
            ey=self.qDif(i,ey)
        return [n,y]

    def KashRep(self,i,w, string=False):
        if isinstance(i,int): i=self.alphabet[i]
        if isinstance(w,Word): x=Wordnomial([w])
        else: x=w
        reps=dict()
        maxn=-1
        while not x.isZero():
            m=self.maxKash(i,x)
            y=m[1]
            n=m[0]
            weight=self.weight(i)
            l=max(abs(self.RS.pair(weight,weight))//2,1)
            sign=(-1)**self.parity[i]
            reps[n]=Monomial(sign**(((n-1)*n)//2),l*((n-1)*n)//2)*y
            #print(x)
            #print(y)
            x=x+Monomial(-1,0)*self.shuffle((self.divpow(self.RV[weight],n,length=l,sign=sign)),reps[n])
        m=str(w)+" = "
        init=True
        for n in reps:
            if init: init=not init
            else: m=m+"+ "
            m=m+str(i)+"^("+str(n)+") * ("+str(reps[n])+" ) "
        if not string: return reps
        else: return m
        
    def kashF(self,i,w, p=1, debug=False):
        if isinstance(i,int): i=self.alphabet[i]
        if isinstance(w,Word): st=w.strWord()
        elif isinstance(w,Wordnomial):
            if w.strWords()==[]: st=""
            else: st=w.strWords()[0]
        nw=Wordnomial([])
        r=self.KashRep(i,w)
        if debug: print("kashF Debug:\t r="+str([(x,str(r[x])) for x in r]))
        weight=self.weight(i)
        l=max(abs(self.RS.pair(weight,weight))//2,1)
        if debug: print("kashF Debug:\t l="+str(l))
        sign=(-1)**self.parity[i]
        for n in r:
            rv=self.RV[weight]
            if debug: print("kashF Debug:\t Root Vector="+str(rv))
            dp=self.divpow(self.RV[weight],n+1,length=l,sign=sign, debug=debug)
            if debug: print("kashF Debug:\t "+str(n+1)+"th div pow="+str(dp))
            nw=nw+self.shuffle(dp,r[n])
        if p==1: return nw
        else: return self.kashF(i,nw,p=p-1)

    def tikzCrystal(self,depth=1):
        rank=self.RS.rank
        tikz=get_N_colors(rank)
        tikz+="\n\\begin{tikzpicture}[scale=1]\n"
        def f(i, x):
            fx=self.PBWexp(self.kashF(i,x))
            fxmodq=fx.modQ()
            if fxmodq.isZero(): return 0
            s=fxmodq.listRep()

            coef=s[0].coefficient
            if len(s)==1: 
                w=s[0]
                coef=w.coefficient()
                if coef==1 : return (1,w.getWord()[1])
                elif coef==-1: return (-1,w.getWord()[1])
                else:
                    raise ArithmeticError("PBW expansion of f_"+str(i)+" "+str(self.PBWexp(x))+" is "+str(fx))
            else:
                raise ArithmeticError("PBW expansion mod q has more than one term:"+str([str(a) for a in s]))


        words={0:[""]}
        tikz+="\\draw (0,0) node (0){$\\emptyset$};\n"

        maxheight=0
        hunit=2.5
        vunit=.5
        signdict={1:"", -1:"diamond"}
        for d in range(1,depth):
            temp=[]
            comp=compositions(rank-1,d)
            for x in [Vector(c) for c in comp]:
                temp+=self.PBWwords(x)
            words[d]=sorted(temp[:])

            count=0
            maxheight=(len(words[d])-1)*vunit
            for w in words[d]:
                tikz+="\\draw "+str((hunit*d,maxheight-(2*vunit)*count))+" node ("+w+") {$"+w+"$};\n"
                count+=1

            for w in words[d-1]:
                for i in range(rank):
                    fw=f(i,self.PBWel(w))
                    if w=="":pw="0"
                    else: pw=w
                    if fw!=0:
                        tikz+="\\draw[color=color"+str(i)+","+signdict[fw[0]]+"->, thick,] ("+str(pw)+".east)  -- ("+str(fw[1])+".west);\n"
                count+=1
        tikz+="\\end{tikzpicture}\n\n"
        print(tikz)

    def gojsCrystal(self,depth=1):
        json='{ "class": "go.GraphLinksModel",'
                
        rank=self.RS.rank
        colors=gojs_get_N_colors(rank)
        nodes='\t"nodeDataArray": [ '
        links='\t"linkDataArray": [ '
        def f(i, x):
            fx=self.PBWexp(self.kashF(i,x))
            fxmodq=fx.modQ()
            if fxmodq.isZero(): return 0
            s=fxmodq.listRep()

            coef=s[0].coefficient
            if len(s)==1: 
                w=s[0]
                coef=w.coefficient()
                if coef==1 : return (1,w.getWord()[1])
                elif coef==-1: return (-1,w.getWord()[1])
                else:
                    raise ArithmeticError("PBW expansion of f_"+str(i)+" "+str(self.PBWexp(x))+" is "+str(fx))
            else:
                raise ArithmeticError("PBW expansion mod q has more than one term:"+str([str(a) for a in s]))


        words={0:[""]}
        nodes+='\t{ "key":"", "text": "" },\n'

        maxheight=0
        hunit=2.5
        vunit=.5
        signdict={1:"", -1:"diamond"}
        for d in range(1,depth):
            temp=[]
            comp=compositions(rank-1,d)
            for x in [Vector(c) for c in comp]:
                temp+=self.PBWwords(x)
            words[d]=sorted(temp[:])

            count=0
            maxheight=(len(words[d])-1)*vunit
            for w in words[d]:
                nodes+='\t{ "key":"'+w+'", "text": "'+w+'" },\n'
                count+=1
            for w in words[d-1]:
                for i in range(rank):
                    fw=f(i,self.PBWel(w))
                    if w=="":pw="0"
                    else: pw=w
                    if fw!=0:
                        if fw[0]==1: sgn=""
                        else: sgn="-"
                        links+='\t{ "from": "'+w+'", "to": "'+fw[1]+'", "color": '+colors[i]+', "text": "'+sgn+str(self.alphabet[i])+'" },\n'
                count+=1
        nodes=nodes[:-2]+"],\n"
        links=links[:-2]+"]\n"
        json+=nodes+links+"}"
        return json

    def genDCB(self,weight):
        if weight in self.DCB:
            return self.DCB[weight]
        dCB=[]
        wordpow=self.PBWwordpowDict(weight)
        words=[]
        for x in wordpow:
            words+=[x]
        words=sorted(words, key=self.lexKey)
        for i in range(len(words)):
            w=words[i]
            b=self.PBWel(w,dual=True)
            #print((w,str(b)))
            if b[w].coefficient().nonBarSym(aniso=False).isZero(): 
                aniso=False
            else: 
                aniso=True
            #print((str(b),aniso))
            if i!=0:
                less=words[i-1::-1]
                for w2 in less:
                    ww2=self.weight(w2)
                    lower=self.DCB[w2]
                    cbcoef=b[w2].coefficient().nonBarSym(aniso=aniso)
                    ocbcoef=cbcoef
                    lowercoef=lower[w2].coefficient()
                    lowerdeg=lowercoef.maxDeg()
                    mult=Polynomial([])
                    count=0
                    while not cbcoef.isZero():
                        md=cbcoef.maxDeg()
                        ratio=cbcoef[md].coefficient()//(lowercoef[lowerdeg].coefficient())
                        if md==lowerdeg:
                            cbcoef=(cbcoef-(ratio//abs(ratio))*Monomial.Q(md-lowerdeg)*lowercoef).nonBarSym(aniso=aniso)
                            mult+=(ratio//abs(ratio))*Monomial.Q(md-lowerdeg)
                        else: 
                            cbcoef=(cbcoef-ratio*Monomial.Q(md-lowerdeg)*lowercoef).nonBarSym(aniso=aniso)
                            mult+=ratio*Monomial.Q(md-lowerdeg)
                        count+=1
                        if md<=0 or count>6000:
                            print("Calculation of multiple of PBW^*("+w2+") in CB^*("+w+") resulted in invalid multiple")
                            print(words) 
                            printN(aniso)
                            printN(count)
                            printN(w+"<-"+w2)
                            printN(str(self.PBWexp(self.PBWel(w,dual=True),dual=True))+"="+str(self.PBWel(w,dual=True)))
                            printN(str(self.PBWexp(self.PBWel(w2,dual=True),dual=True))+"="+str(self.PBWel(w2,dual=True)))
                            printN(lower)
                            printN(b[w2].coefficient())
                            printN(b[w2].coefficient().nonBarSym())
                            printN(b[w2].coefficient().nonBarSym(aniso))
                            printN(b[w2].coefficient().nonBarSym(aniso=aniso))
                            printN(ocbcoef)
                            printN(cbcoef)
                            printN(lowercoef)
                            printN(cbcoef-ratio*Monomial.Q(md-lowerdeg)*lowercoef)
                            printN((cbcoef-ratio*Monomial.Q(md-lowerdeg)*lowercoef).nonBarSym(aniso))
                            quit()
                    b=b-mult*lower
            #print("DCB["+w+"]="+str(b) +" = "+str(self.PBWexp(b, dual=True)))
            dCB+=[b]
            self.DCB[w]=b

        self.DCB[weight]=dCB
        return dCB

    def genCB(self,weight):
        if weight in self.CB:
            return self.CB[weight]
        dCB=self.genDCB(weight)
        dCBexp={}
        CB=[]
        wordpow=self.PBWwordpowDict(weight)
        words=[]
        for x in wordpow:
            words+=[x]
            dCBexp[x]={}
            dexp=self.PBWexp(self.DCB[x], dual=True)
            for i in dexp.listRep():
                dCBexp[x][i.getWord()[1]]=i.coefficient()
            #print([(y, str(dCBexp[x][y])) for y in dCBexp[x]])
        words=sorted(words, key=self.lexKey)
        numwords=len(words)
        for i in range(numwords):
            coefs={}
            w1=words[i]
            cb=Wordnomial([])
            for j in range(numwords):
                w2=words[j]
                if i==j: pair=1
                else: pair=0
                coefs[j]=Polynomial([Monomial(pair,0)])
                for k in range(j):
                    w3=words[k]
                    if w3 in dCBexp[w2]: coefs[j]=coefs[j]-dCBexp[w2][w3]*coefs[k]
                cb+=coefs[j]*self.PBWel(w2)
            if self.bar(cb)!=cb: 
                raise ValueError("Error: Canonical basis element is not bar invariant:\n"+ 
                    "CB["+w1+"]="+str(self.PBWexp(cb))+"\n"+
                    str(cb)+"\n"+str(self.bar(cb)))
            CB+=[cb]
            self.CB[w1]=cb
        self.CB[weight]=CB
        return CB


    def CBexp(self,x, weight=None, string=False, dual=False, debug=False): #x is a wordnomial of weight a alpha_1 + b alpha_2 + c alpha_3
        if debug: print("CBexp in debug mode")
        if weight is None:
            wd=self.weightDecomp(x)
            weights=wd.keys()
            s=""
            expansion=Wordnomial([])
            initial=True
            for weight in weights:
                if not initial: 
                    s+= " + "
                else:
                    initial=False
                if string: s+=self.CBexp(wd[weight],weight, string, dual, debug)
                else:
                    cbe=self.CBexp(wd[weight],weight, string, dual, debug)

                    try:
                        expansion+=cbe
                    except TypeError:
                        raise TypeError("Expansion"+str(cbe))
            if string: return s
            else: return expansion
        self.genCB(weight)
        domW=self.PBWwords(weight)
        CBmult=dict()
        y=self.PBWexp(x,weight,False,dual,debug)
        domW=sorted(domW,key=self.lexKey)
        if debug: print(domW)
        m=""
        expansion=Wordnomial([])
        prefix="PBW"
        if dual: prefix+="^*"
        prefix+=":"
        if dual:
            CB=self.DCB
        else:
            CB=self.CB
        for w in domW:
            b=self.PBWexp(CB[w],weight,False,dual,debug)
            if isinstance(y,str):
                #print(y)
                y=Wordnomial([Word(Monomial.ONE(),y)])
            if isinstance(y,Word): y=Wordnomial([y])
            wcoef=(y[prefix+w]).coefficient()
            cbcoef=b[prefix+w].coefficient()
            if cbcoef.isZero(): 
                raise ArithmeticError("Coefficient of PBW("+w+") in CB("+w+") is zero!?")
            CBmult[w]=(wcoef/cbcoef)[0]
            y=y-(CBmult[w]*b)

        init=True
        for w in domW:
            prefix="CB"
            if dual: prefix+="^*"
            expansion=expansion+Word(CBmult[w],[prefix+":",w])
            if not CBmult[w].isZero():
                if init: init=not init
                else: m=m+"+ "
                if (CBmult[w]-Monomial(1,0)).isZero(): m=m+prefix+"("+w+")"
                else: m=m+"( "+str(CBmult[w])+" ) "+prefix+"("+w+")"
        if y.isZero(): 
            if string: return m
            else: return expansion
        else: 
            print("Error: "+str(x)+"=("+m+")+"+str(y))
            quit()



U=QuantumGroup(RootSystem.buildByType("G", [0,0], 1), useCache=True, debug=False)
shf=U.shuffle
pbw=U.PBWel
mon=U.MonEl
f=U.kashF
print(U.LynStr)
#U=QuantumGroup(RootSystem.buildByType("F", [0,0,0,0], 1), useCache=True, debug=True)
#print(U.rootLattice(Vector([2,-2])))
#U.gojsCrystal(depth=3)
"""for x in U.RV:
    print((str(U.RVmon[x][0]),str(U.RV[x])))

for d in range(1,8):
    comp=compositions(U.rank-1,d)
    for v in [Vector(c) for c in comp]:
        U.genCB(v)
        words=U.PBWwords(v)
        print("Finished generating CB for "+str(v))
        for w in words:
            print(("CB["+w+"]="+str(U.PBWexp(U.CB[w]))).ljust(100)+("PBW["+w+"]="+str(U.CBexp(U.PBWel(w)))).ljust(50))"""
#U.appendCache(debug=True)