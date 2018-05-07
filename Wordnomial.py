##############################################################
# Wordnomial
# Author: Sean Clark
#
# A collection of classes and methods to compute
# shuffle products associated to root data.
##############################################################
#!/usr/bin/env python
import math


class Monomial:
    """Basic unit of an integer polynomial; consists of a power of an indeterminant (q) and an integer coefficient."""

     
    def __init__(self,coef,p):   
        #Initializes the integer coefficient and the integer power
        if not isinstance(coef,int) or not isinstance(p,int):
            raise TypeError("Monomial must be initialized with two integer arguments: (coefficient,power); you entered ("
            + str(coef)+", "+str(p)+" )")
        self.power=p
        self.scalar=coef
        self.simplify()

    @classmethod
    def ONE(cls):
        return Monomial(1,0)
    @classmethod
    def ZERO(cls):
        return Monomial(0,0)
    @classmethod
    def Q(cls, n=1):
        return Monomial(1,n)
    def QI(cls):
        return Monomial(1,-1)

    def deg(self):    
        #returns the power
        return self.power
    def bar(self):    
        #returns the power
        return Monomial(self.scalar,-self.power)

    def coefficient(self):     
        #Returns the coefficient
        return self.scalar

    def simplify(self): 
        #Replaces 0*q^n with 0*q^0, and returns self.
        if self.scalar== 0: self.p=0
        return self

    def __mul__(self,other):    
        #Defines multiplication operation with integers and other monomials
        
        if isinstance(other,int): 
            #if integer, scale coefficient and leave power the same.
            newp=self.power
            newcoef=other*self.scalar
            return Monomial(newcoef,newp)
        
        elif isinstance(other, Monomial): 
            #if monomial, multiply coefficients and add powers.
            p=other.deg()
            s=other.coefficient()
            newp=self.power+p
            newcoef=s*self.scalar
            return Monomial(newcoef,newp)
        else:
            return NotImplemented

    def __rmul__(self,other):    
        #Defines right multiplication operation with integers and other monomials
        
        return self*other


    def __add__(self,other): 
        #Defines addition operation with integers, monomials
        
        if isinstance(other,int): 
            #If an integer, either add to coefficient or make polynomial
            if self.power==0:
                newcoef=self.scalar+other
                return Monomial(newcoef,0)
            else:
                mon=Monomial(other,0)
                return Polynomial([self,mon])
        elif isinstance(other,Monomial): 
            #If another monomial, either add coefficients or make polynomial
            p=other.deg()
            s=other.coefficient()
            if self.power==p or self.scalar==0 or s==0:
                if self.scalar==0: newp=p
                else: newp=self.power
                newcoef=s+self.scalar
                return Monomial(newcoef,newp)
            else:
                return Polynomial([self, other])
        

        else: 
            return NotImplemented
    
    def __radd__(self,other): 
        return self+other


    def __sub__(self,other): 
        #Defines subtraction operation with: integers, monomials, polynomials, words, and wordnomials
        
        if isinstance(other,int) or isinstance(other,Monomial):
            return self + (other * (-1))

        else:
            return NotImplemented

    def __rsub__(self,other): 
        #Defines subtraction operation with: integers, monomials, polynomials, words, and wordnomials
        return (self*(-1))+other

    def __str__(self): 
        #Defines a String representation of the form c*q^p
        s=""
        if self.scalar==0:
            return "0"
        elif self.scalar==1 and self.power==0: return "1"
        elif self.scalar==1: s=s+""
        elif self.scalar==-1 and self.power==0: return "-1"
        elif self.scalar==-1: s=s+"-"
        else: s=s+str(self.scalar)
        if self.power==0: s=s+""
        else: s=s+"q^{"+str(self.power)+"}"
        return s

    def __eq__(self, other): 
        #Defines shallow equality: monomials are equal if their coefficients and powers are the same
        #Also allows to check against Polynomial and Wordnomial representations of a monomial
        if isinstance(other,int): 
            return self.power==0 and self.scalar==other
        elif isinstance(other,Monomial):
            return self.scalar==other.coefficient() and self.power==other.deg()
        elif isinstance(other, Polynomial):
            return other.isMonomial() and self==other.listRep()[0]
        elif isinstance(other, Word):
            return len(other)==0 and other.coefficient()==self
        elif isinstance(other, Wordnomial):
            return other.isWord() and self==other.listRep()[0]
        else:
            return False
            
    def cacheStr(self):
        if self.scalar==0:
            return "0"
        s=str(self.scalar)+"*q^"+str(self.power)
        return s
    @classmethod
    def uncache(self,s):
        if s=="0": return Monomial(0,0)
        scal,power=s.split("*")
        return Monomial(int(scal),int(power[2:]))

    @classmethod
    def testSuite(self):
        print("Initialization test:\n--------------------\n")
        ma=Monomial(2,4)
        mb=Monomial(1,3)
        mc=Monomial(-1,-3)
        md=Monomial(0,2)
        try:
            me=Monomial(3,'a')
        except TypeError:
            print("Invalid initialization thrown as intended for Monomial(3,'a')")
        print("Monomials initialized successfully!")

        print("\nString method tests:\n--------------------\n")
        print("Does str(Monomial(2,4)) yield the string 2q^{4}? "+ str((str(ma)=="2q^{4}")))
        print("Does str(Monomial(1,3)) yield the string q^{3}? "+ str((str(mb)=="q^{3}")))
        print("Does str(Monomial(-1,-3)) yield the string -q^{3}? "+ str((str(mc)=="-q^{-3}")))
        print("Does str(Monomial(0,-2)) yield the string 0? "+ str((str(md)=="0")))

        print("\nData method tests:\n--------------------\n")
        print("Does Monomial(2,4).deg() yield 4? "+ str(ma.deg()==4))
        print("Does Monomial(1,3).coefficient() yield 1? "+ str(mb.coefficient()==1))

        print("\nMultiplication tests:\n--------------------\n")
        mbc=mb*mc
        tmc=mc*2
        print("Does Monomial(1,3)*Monomial(-1,-3) have degree 0 and coefficient -1? "
            + str((mbc.deg(),mbc.coefficient())==(0,-1)))
        print("Does Monomial(1,3)*Monomial(-1,-3) leave the factors unchanged? "+
            str(mb.deg()==3 and mb.coefficient()==1 and mc.deg()==-3 and mc.coefficient()==-1))
        print("Does 2*Monomial(-1,-3) have degree -3 and coefficient -2? "+
            str(tmc.deg()==-3 and tmc.coefficient()==-2))
        print("Does 2*Monomial(-1,-3) leave Monomial(-1,-3) unchanged? "+
            str(mc.deg()==-3 and mc.coefficient()==-1))
        try:
            badm=ma*[]
        except TypeError:
            print("Successfully caught invalid multplication of Monomial and List")

#Monomial.testSuite()










class Polynomial:
    #Representation of a sum of monomials using a list of monomials

    def __init__(self,mons):
        #Initializes a polynomial; the argument should be a list of monomials, though some
        #leeway is built in. The init automatically simplifies the polynomial to avoid
        #repeated calls in the methods, and then generates some statistics.
        if isinstance(mons,Monomial):
            #If the argument is monomial, converts to a list of one monomial
            b=[mons]
            mons=b
        if isinstance(mons,int):
            #If the argument is integer, converts to a list of one monomial
            b=[Monomial(mons,0)]
            mons=b
        if not isinstance(mons,list): 
            raise TypeError("Polynomial instantiation error: argument must be a list of monomials")
        valid=True
        for x in mons:
            valid=valid and isinstance(x,Monomial)
        if valid: 
            self.mons=mons
        else:
            raise TypeError("Polynomial instantiation error: argument must be a list of monomials")
        s=[]
        for x in self.mons:
            d=x.deg()
            if not d in s: s=s+[d]
        s.sort()
        self.degrees=s[:]
        self.simplify()
        d=dict()
        for x in self.mons:
            d[x.deg()]=x
        self.degrees=list(d.keys())
        self.degrees.sort()
        self.dict=d

    @classmethod
    def ONE(cls):
        return Polynomial(Monomial.ONE())
    @classmethod
    def ZERO(cls):
        return Polynomial(Monomial.ZERO())

    @classmethod
    def qInt(cls,n, sign=1, length=1, debug=False):#Produces a quantum integer
        ml=[]
        if length==0: return Polynomial.ONE()
        if debug: 
            print("qInt debug: n="+str(n)+"\t sign="+str(sign)+"\t length="+str(length))
        for i in range(0,n):
            ml=ml+[Monomial(sign**(n-1-i),(n-1-2*i)*length)]
        return Polynomial(ml)

    @classmethod
    def qFactorial(cls,n, sign=1, length=1):#Produces a quantum integer
        if not (isinstance(n, int) and n>=0): raise TypeError('Quantum factorial only defined for non-negative integer n')
        if n==0 or n==1: return Polynomial.qInt(1, sign, length)
        else:
            return Polynomial.qInt(n, sign, length)*Polynomial.qFactorial(n-1, sign, length)

    @classmethod
    def qBin(cls,n, a, sign=1, length=1):#Produces a quantum binomial
        num=Polynomial.qInt(1)
        den=Polynomial.qInt(1)
        for i in range(a):
            num*=(Monomial(sign**(n-i),(n-i)*length)-Monomial(1,(i-n)*length))
            den*=(Monomial(sign**(1+i),(1+i)*length)-Monomial(1,(-1-i)*length))
        return (num/den)[0]

    def simplify(self):
        #Method for combining like terms to simplify the polynomial, and removing zero terms.
        smons=[]
        if self.mons==[]: return self
        else:
            for d in self.degrees:
                nm=Monomial(0,d)
                for m in self.mons:
                    if m.deg()==d: nm=nm+m
                if(not nm.coefficient()==0): smons=smons+[nm]
        self.mons=smons
        return self

    def __getitem__(self,k):
        #Finds the monomial of degree k if it exists; else returns 0
        try:
            if isinstance(k,int): return self.dict[k]
            elif isinstance(k,slice):
                start=k.start
                if start is None: start=self.minDeg()
                stop=k.stop
                if stop is None: stop=self.maxDeg()+1
                step=k.step
                if step is None: step=1
                l=[]
                for i in range(start,stop,step):
                    if i in self.dict:
                        l+=[self.dict[i]]
                return Polynomial(l)
            else: raise KeyError()

        except KeyError:
            return Monomial(0,0)

    
    def listRep(self):
        #returns a duplicate of the list representation of the polynomial
        return self.mons[:]

    
    def dictRep(self):
        #returns a duplicate of the dict representation of the polynomial: (degree:monomial)
        return self.dict.copy()

    def __len__(self):
        #Returns the number of summands 
        return len(self.mons)

    def maxDeg(self):
        #finds the largest degree of a monomial with nonzero coefficient
        if len(self.degrees)==0: return 0
        else: return self.degrees[-1]

    def minDeg(self):
        #finds the smallest degree of a monomial with nonzero coefficient
        if len(self.degrees)==0: return 0
        else: return self.degrees[0]

    def leadingTerm(self):
        #returns the monomial of highest degree
        if len(self.degrees)==0: return Monomial(0,0)
        else: return self.dict[self.degrees[-1]]

    def barSym(self,length=1):
        sym=Polynomial([])
        if length>0:
            m=max(-self.minDeg(),self.maxDeg())
            for d in range(0,m+1):
                if d==m:
                    t1=(self[-d].bar())
                    t2=self[d]
                sym+=(self[-d].bar()-self[d])
        return sym

    def nonBarSym(self,length=1, aniso=False):
        sym=Polynomial([])
        if length>0:
            if aniso: sign=-1
            else: sign=1
            m=max(-self.minDeg(),self.maxDeg())
            for d in range(0,m+1):
                if d==m:
                    t1=(self[-d].bar())
                    t2=self[d]
                sym+=(self[d]-(sign**d)*self[-d].bar())
        return sym

    def bar(self):
        new=Polynomial([])
        for d in self.degrees:
            new+=self.dict[d].bar()
        return new


    def deg(self):
        #A shorter command to be consistent with Monomial.
        return self.maxDeg()

    def degs(self):
        #Returns the set of degrees of the nonzero monomials in increasing order
        return self.degrees[:]

    def isZero(self):
        #Checks for equality with zero, represented in this class as an empty sum;
        return self.mons==[]

    def isMonomial(self):
        #Checks if the polynomial is simply one monomial
        self.simplify()
        return len(self.mons)==1
        
    def __add__(self, other):
        #Defines addition for polynomials

        if isinstance(other,int):
            #if an integer, convert to a monomial and add to the list
            newmons=self.mons+[Monomial(other,0)]
            poly= Polynomial(newmons)
            return poly.simplify()
        elif isinstance(other,Monomial):
            #if a monomial, add to the list
            newmons=self.mons+[other]
            poly= Polynomial(newmons)
            return poly.simplify()
        elif isinstance(other,Polynomial):
            #if a polynomial, concatenate their lists
            newmons=self.mons+other.listRep()
            poly= Polynomial(newmons)
            return poly.simplify()
        elif isinstance(other,Word) or isinstance(other, Wordnomial):
            #if a word, convert the polynomial to a Wordnomial and use Wordnomial addition
            new=Wordnomial([Word(self,[])])
            return new+other
        else:
            raise TypeError("Polynomials may only be added to wordnomials, words, polynomials, monomials, or integers.")


    def __sub__(self,other): 
        #Defines subtraction operation with: integers, monomials, polynomials, words, and wordnomials
        
        if (isinstance(other,int) or isinstance(other,Monomial) or isinstance(other, Polynomial) 
            or isinstance(other, Word) or isinstance(other, Wordnomial)):
            return self + (other * -1)

        else:
            raise TypeError("Monomials may only be added to Integers, Polynomials, Words, or Wordnomials")

    def __mul__(self, other):
        mons=self.mons
        newmons=[]
        if isinstance(other,int) or isinstance(other,Monomial):
            #If an integer, scale all the summands
            for x in mons:
                newmons=newmons+[x*other]
            poly= Polynomial(newmons)
            return poly

        elif isinstance(other,Polynomial):
            #if a polynomial, multiply entry by entry
            otherm=other.listRep()
            for om in otherm:
                for m in mons:
                    nm=om*m
                    newmons=newmons+[nm]
            poly= Polynomial(newmons)
            return poly

        elif isinstance(other,Word) or isinstance(other, Wordnomial):
            #if a Word or Wordnomial, use that multiplication
            return other*self

        else:
            raise TypeError("Polynomials may only multiply wordnomials, words, polynomials, monomials, or integers.")
    def __rmul__(self, other):
        return self*other
            
    def __truediv__(self, other): # polynomial division, returns a pair (quotient, remainder)
        if isinstance(other,int) or isinstance(other,Monomial):
            temp=Polynomial(other)
            other=temp
        if other.isZero():
            raise ArithmeticError("Division by zero!")
        nself=Monomial(1,-self.minDeg())*self
        nother=Monomial(1,-other.minDeg())*other
        fact=Monomial(1,-other.minDeg()+self.minDeg())
        if nself.maxDeg()<nother.maxDeg():
            return (Polynomial(0), self)
        
        quo=Polynomial(0)
        rem=nself
        lcdiv=nother.leadingTerm().coefficient()
        degdiv=nother.maxDeg()
        degn=nself.maxDeg()
        for p in range(0,degn+1-degdiv):
            c=rem[degn-p].coefficient()
            if (c==0) or (c%lcdiv!=0):
                mul=Polynomial(0)
            else:
                mul=Monomial(c//lcdiv,degn-degdiv-p)
            quo=quo+mul
            rem=rem-mul*nother
            #print(str(p)+"  "+str(degdiv)+"  "+str(c)+"  "+str(mul)+"  "+str(quo)+"  "+str(rem))
        return (fact*quo,rem)

    def __div__(self, other): # polynomial division, returns a pair (quotient, remainder)
        return self.__truediv__(other)

    def __str__(self):
        #Provides a string representation of the polynomial as a sum of monomials
        s=""
        if self.isZero():
            s="0"
        else:
            init=True #When init is true, no + sign is prefixed onto the monomial string
            for m in self.mons:
                if init:
                    s=s+str(m)
                    init=False
                else: 
                    if m.scalar>0:

                        s=s+" + " + str(m)
                    else:
                        s=s+" - "+str(-1*m)
        return s
            
    def cacheStr(self):
        if self.isZero():
            return "0"
        else:
            s=""
            init=True #When init is true, no + sign is prefixed onto the monomial string
            for m in self.mons:
                if init:
                    s+="("+m.cacheStr()+")"
                    init=False
                else: 
                    s+="+("+m.cacheStr()+")"

        return s
    @classmethod
    def uncache(self,s):
        if s=="0": return Polynomial([])
        l=s.split("+")
        ml=[]
        for x in l:
            ml+=[Monomial.uncache(x[1:-1])]
        return Polynomial(ml)
        
    def __eq__(self, other): 
        #Defines shallow equality
        if isinstance(other,int): 
            return self.isMonomial() and self.mons[0]==other
        elif isinstance(other,Monomial):
            return self.isMonomial() and self.mons[0]==other
        elif isinstance(other, Polynomial):
            degs=self.degrees
            if degs==other.degs():
                if self.isZero(): return other.isZero()
                equal=True
                for x in degs:
                    equal=equal and self[x]==other[x]
                return equal
            else: return False
        elif isinstance(other, Word):
            return len(other)==0 and other.coefficient()==self
        elif isinstance(other, Wordnomial):
            return other.isWord() and self==other.listRep()[0]
        else:
            return False

    @classmethod
    def testSuite(self):
        pa=Polynomial([Monomial(4,2),Monomial(-2,2),Monomial(-1,0),Monomial(2,0)])
        pb=Polynomial([Monomial(2,2), Monomial(1,0)])
        print("Does 4q^2-2q^2-1+2=2q^2+1?"+
            str(pa==pb))
        print("Does 4q^2-2q^2-1+2 have leading term 2q^2? "+
            str(pa.leadingTerm()==Monomial(2,2)))
        print("Does 4q^2-2q^2-1+2-2q^2 have leading term 1? "+
            str(pa-Monomial(2,2).simplify()==Monomial(1,0)))
        pc=Polynomial([Monomial(1,1),Monomial(-1,0)])
        pd=Polynomial([Monomial(2,2), Monomial(1,0)])
        x=pc/pd
        print(str(pc)+" = ("+ str(x[0])+") * ( "+str(pd)+" ) + ( "+str(x[1])+" )")
        pc=Polynomial([Monomial(1,2),Monomial(-1,-2)])
        pd=Polynomial([Monomial(1,1), Monomial(-1,-1)])
        x=pc/pd
        print(str(pc)+" = ("+ str(x[0])+") * ( "+str(pd)+" ) + ( "+str(x[1])+" )")





class Word:
    """Akin to Monomial, represents a word in some alphabet along with some polynomial coefficient."""

    def __init__(self,coef,w):
        if isinstance(coef,int): coef=Polynomial([Monomial(coef,0)])
        if isinstance(coef,Monomial): 
            coef=Polynomial([coef])
        if isinstance(w,str): w=w.split('')
        if not isinstance(w,list) or not isinstance(coef,Polynomial):
            raise TypeError("Word takes arguments of the form (Polynomial, List)")
        self.word=w
        self.scalar=coef
        self.simplify()

    def simplify(self):
        if self.scalar.isZero(): self.word=[]
        return self

    def isZero(self):
        return self.scalar.isZero()



    @classmethod
    def Empty(cls):
        return Word(Monomial(1,0),[])

    def divQ(self, length=1):
        if length>0: return self.scalar.minDeg()>=length
        elif length<0: return self.scalar.maxDeg()<=length

    def modQ(self, length=1):
        if length>0: return Word(self.scalar[:length],self.word)
        elif length<0: return Word(self.scalar[length+1:],self.word)

    def __getitem__(self,w):
        #Returns the coefficient of the given word w if it exists; else returns 0
        if w==self.strWord():
            return self
        else:
            return Word(Polynomial([Monomial(0,0)]),[])



    def getWord(self):
        #Returns the list representation of the word
        return self.word[:]

    def revWord(self):
        #Reverses the word, essentially applying an anti-involution
        return Word(self.scalar,self.word[::-1])

    def strWord(self):
        #Returns the string representation of the word
        s=""
        for x in self.word: s=s+str(x)
        return s

    def isEmpty(self):
        #Declares if the word is a scalar, i.e. a multiple of the empty word
        return len(self)==0

    def compare(self,other,lex=None):
        if isinstance(other, Word):
            if lex==None:
                if self.getWord()<other.getWord(): return -1
                elif self.getWord()==other.getWord(): return 0
                else: return 1
            oLen=len(other)
            sLen=len(self)
            l=min(oLen,sLen)
            oWord=other.getWord()
            sWord=self.word
            for i in range(l):
                sVal=lex[sWord[i]]
                oVal=lex[oWord[i]]
                if sVal<oVal: return -1
                if sVal>oVal: return 1
            if sLen<oLen: return -1
            elif sLen==oLen: return 0
            else: return 1
        else:
            return NotImplemented

    def __len__(self):
        #Returns the number of letters in the word
        return len(self.word)

    def coefficient(self):
        #Returns the coefficient of the word
        return self.scalar

    def __mul__(self,other):
        #Defines multiplication of words with wordnomials, words, polynomials, etc.

        if isinstance(other,int) or isinstance(other,Monomial) or isinstance(other,Polynomial):
            #If an integer, monomial, or polynomial, then scale the coefficient 
            c=self.scalar*other
            w=self.word[:]
            return Word(c,w)

        elif isinstance(other, Word):
            #If another word, multiply the coefficients and concatenate the words
            #print("Multiplying words")
            c=other.coefficient()*self.scalar
            #print(other.coefficient())
            w=self.word+other.getWord()
            return Word(c,w)
        
        elif isinstance(other, Wordnomial):
            #If a wordnomial, defer to wordnomial multiplication
            s=Wordnomial([self])
            return s*other
        else:
            return NotImplemented

    def __rmul__(self, other):

        if isinstance(other,int) or isinstance(other,Monomial) or isinstance(other,Polynomial):
            return self*other
        if isinstance(other, Word):
            #If another word, multiply the coefficients and concatenate the words
            #print("Multiplying words")
            c=other.coefficient()*self.scalar
            #print(other.coefficient())
            w=other.getWord()+self.word
            return Word(c,w)
        else: 
            return NotImplemented

    def __add__(self,other):
        #Defines addition of words with wordnomials, words, polynomials, etc.
        if isinstance(other,int):
            w=self.word[:]
            new= Word(Polynomial(Monomial(other,0)),w)
            return new
        elif isinstance(other,Monomial):
            w=self.word[:]
            new= Word(Polynomial(other),w)
            return new
        elif isinstance(other,Polynomial):
            w=self.word[:]
            new= Word(other,w)
            return new
        elif isinstance(other, Word):
            new=Wordnomial([self,other])
            return new
        elif isinstance(other, Wordnomial):
            new=Wordnomial([self])+other
            return new
        else:
            raise TypeError("Words may only be added with wordnomials, polynomials, etc.")



    def __sub__(self,other): 
        #Defines subtraction operation with: integers, monomials, polynomials, words, and wordnomials
        
        if (isinstance(other,int) or isinstance(other,Monomial) or isinstance(other, Polynomial) 
            or isinstance(other, Word) or isinstance(other, Wordnomial)):
            return self + (other * -1)

        else:
            raise TypeError("Monomials may only be added to Integers, Polynomials, Words, or Wordnomials")

    def __eq__(self,other):
        return (self-other).isZero()
    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        #Defines a string representation of the word
        s=""
        if self.scalar.isZero():   
            return "0"

        else:
            l=len(self.scalar)
            if l>1: s="("+str(self.scalar)+") "
            else:
                coef=self.scalar.listRep()[0]
                if coef==1: s=s+""
                elif coef==-1: s=s+"-"
                elif coef==0: 
                    return "0"
                else: s=str(coef)+" "
            s=s+ "["
            for i in self.word: s=s+i
            s=s+"]"
            return s
            
    def cacheStr(self):
        if self.isZero():
            return "0"
        else:
            return "("+self.scalar.cacheStr()+")["+self.strWord()+"]"
    @classmethod
    def uncache(self,s):
        if s=="0": return Word(Polynomial([]),[])
        poly,word=s.split(")[")
        return Word(Polynomial.uncache(poly[1:]),list(word[:-1]))




class Wordnomial:

    def __init__(self, s):
        #Initializes a wordnomial; the argument should be a list of words, though some
        #leeway is built in. The init automatically simplifies the polynomial to avoid
        #repeated calls in the methods, and then generates some statistics.

        if isinstance(s,str):
            s.split("+")

        if isinstance(s,Word):
            #If the argument is monomial, converts to a list of one word
            b=[s]
            s=b
        if isinstance(s,Polynomial):
            #If the argument is monomial, converts to a list of one word
            b=[Word(s,[])]
            s=b
        if isinstance(s,Monomial):
            #If the argument is monomial, converts to a list of one word
            b=[Word(Polynomial([s]),[])]
            s=b
        if isinstance(s,int):
            #If the argument is integer, converts to a list of one word
            b=[Word(Polynomial([Monomial(mons,0)]),[])]
            s=b
        if not isinstance(s,list): 
            raise TypeError("Wordnomial instantiation error: argument must be a list of words")
        valid=True
        for x in s:
            valid=valid and isinstance(x,Word)
        if valid: 
            self.wordnom=s
        else:
            raise TypeError("Wordnomial instantiation error: argument must be a list of words")
        s=[]
        w=[]
        for x in self.wordnom:
            wstring=x.strWord()
            wlist=x.getWord()
            if not wstring in s: s=s+[wstring]
            if not wlist in w: w=w+[wlist]
        self.stringwords=s
        self.listwords=w
        self.simplify()
        d=dict()
        for x in self.wordnom:
            d[x.strWord()]=x
        self.dict=d

    def simplify(self):
        smons=[]
        words=self.getWords()
        for w in words:
            nc=Polynomial([Monomial(0,0)])
            for summand in self.wordnom:
                if summand.getWord()==w:
                    nc=nc+summand.coefficient()
            if not nc.isZero(): smons=smons+[Word(nc,w)]
        self.wordnom=smons
        s=[]
        w=[]
        for x in self.wordnom:
            wstring=x.strWord()
            wlist=x.getWord()
            if not wstring in s: s=s+[wstring]
            if not wlist in w: w=w+[wlist]
        self.stringwords=s
        self.listwords=w
        d=dict()
        for x in self.wordnom:
            d[x.strWord()]=x
        self.dict=d
        return self

    def __getitem__(self,w):
        #Returns the coefficient of the given word w if it exists; else returns 0
        try:
            return self.dict[w]
        except KeyError:
            return Word(Polynomial([Monomial(0,0)]),[])

    def __len__(self,w):
        #Returns the number of words
        return len(self.wordnom)

    def isZero(self):
        if len(self.wordnom)==0: return True
        if self.wordnom[0].coefficient().isZero(): return True
        return False

    def divQ(self, length=1):
        for w in self.wordnom:
            if not w.divQ(length): return False
        return True

    def modQ(self, length=1):
        l=[]
        for w in self.wordnom:
            l+=[w.modQ(length)]
        return Wordnomial(l)

    def listRep(self):
        return self.wordnom[:]

    def getWords(self):
        return self.listwords[:]

    def revWord(self):
        new=Wordnomial([])
        for x in self.wordnom:
            new=new+x.revWord()
        return new

    def strWords(self):
        return self.stringwords[:]

    def isWord(self):
        return len(self.wordnom)==1

    def normalize(self,p): #divides each coefficient by p if possible
        if not isinstance(p,Polynomial): raise TypeError("Must be Polynomial normalization")
        nw=Wordnomial(Word(Polynomial([Monomial(0,0)]),[]))
        for w in self.wordnom:
            word=w.getWord()
            coef=w.coefficient()
            div=coef/p
            if not div[1].isZero(): raise ArithmeticError("Word "+str(w)+" has coefficient not divisible by "+str(p))
            newword=Word(div[0],word)
            nw=nw+newword
        return nw
        

    def __add__(self,other):
        valid=False
        new=self.wordnom[:]
        if isinstance(other,int) or isinstance(other,Monomial) or isinstance(other,Polynomial):
            new= new+[Word(other,[])]
            valid=True
        if isinstance(other, Word):
            new= new+[other]
            valid=True
        if isinstance(other, Wordnomial):
            new=new+other.listRep()
            valid=True
        if not valid:
            return NotImplemented
        return Wordnomial(new)

    def __mul__(self,other):
        valid=False
        new=[]
        if isinstance(other,int) or isinstance(other,Monomial) or isinstance(other,Polynomial):
            if isinstance(other, int): nother=Monomial(other,0)
            else: nother=other
            for each in self.wordnom:
                new=new+[nother*each]
            valid=True
        if isinstance(other, Word):
            for each in self.wordnom:
                new= new+[each*other]
            valid=True
        if isinstance(other, Wordnomial):
            for each in self.wordnom:
                for eachother in other.listRep():
                    new=new+[each*eachother]
            valid=True
        if not valid:
            print("Error: Unable to multiply Word and " + str(type(other)))
            new=Wordnomial([Word(0,[])])
        return Wordnomial(new)

    def __rmul__(self,other):
        valid=False
        new=[]
        if isinstance(other,int) or isinstance(other,Monomial) or isinstance(other,Polynomial):
            return self*other
        if isinstance(other, Word):
            for each in self.wordnom:
                new= new+[other*each]
            valid=True
        if isinstance(other, Wordnomial):
            for each in self.wordnom:
                for eachother in other.listRep():
                    new=new+[eachother*each]
            valid=True
        if not valid:
            print("Error: Unable to multiply Word and " + str(type(other)))
            new=Wordnomial([Word(0,[])])
        return Wordnomial(new)


    def __sub__(self,other): 
        #Defines subtraction operation with: integers, monomials, polynomials, words, and wordnomials
        
        if (isinstance(other,int) or isinstance(other,Monomial) or isinstance(other, Polynomial) 
            or isinstance(other, Word) or isinstance(other, Wordnomial)):
            return self + (other * -1)

        else:
            raise TypeError("Monomials may only be added to Integers, Polynomials, Words, or Wordnomials")

    def __eq__(self,other):
        return (self-other).isZero()

    def __ne__(self,other):
        return not (self==other)
    def __hash__(self):
        return hash(str(self))
            
    def __str__(self):
        if self.wordnom==[]: return "0"
        s=""
        init=True
        for w in self.wordnom:
            if init:
                s=s+str(w)
                init=False
            else:
                s=s+" + "+str(w)
        return s
            
    def cacheStr(self):
        if self.isZero():
            return "0"
        else:
            s=""
            init=True #When init is true, no + sign is prefixed onto the monomial string
            for w in self.wordnom:
                if init:
                    s+="("+w.cacheStr()+")"
                    init=False
                else: 
                    s+=" + ("+w.cacheStr()+")"

        return s
    @classmethod
    def uncache(self,s):
        if s=="0": return Wordnomial([])
        l=s.split(" + ")
        ml=[]
        for x in l:
            ml+=[Word.uncache(x[1:-1])]
        return Wordnomial(ml)
            
