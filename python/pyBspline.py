#!/usr/bin/env python
# coding: utf-8

# In[24]:


get_ipython().system('jupyter nbconvert --to script pyBspline.ipynb')


# ## Shape class

# In[15]:


class shape :
    
    ###
    def __init__(self,dim=0,codim=0):
        self._dim   = dim
        self._codim = codim
        
    ###
    def dim(self)   : return self._dim
    def codim(self) : return self._codim
    
    ###
    def show(self):
        print("dimension   : ",self._dim)
        print("codimension : ",self._codim)
        
    ###
    #def astuple(self):
    #    return (self.dim,self.codim)
        


# ## Knot vector class

# In[16]:


import numpy as np

class knot_vector:
    
    #_pol_order  = 0 #grado polinomiale
    #_basis_card = 0 #cardinalità della base 
    #_vect       = np.ndarray(shape=(0,))
      
    ###
    def __init__(self,
                 p=0,
                 n=0,
                 cls=np.zeros(shape=(0,1)),
                 check_open=True):
        self._vect       = np.asarray(cls)
        self._pol_order  = p
        self._basis_card = n
        
        if p+n+1 != len(cls) :
            print ("knot_vector")
            print ("p : " , p )
            print ("n : " , n )
            print ("v : " , len(cls) )
            print ("knot_vector error : wrong vector size")
            raise Exception()
            
        # checking knot vector is sorted
        if not np.all(np.diff(self._vect) >= 0) :
            raise Exception("knot vector not sorted")
            
        if check_open == True and self.is_open() == False :
            print ("knot_vector error : it is not an OPEN knot vector")
            raise Exception()
    
    ###
    def p     (self): return self._pol_order
    def n     (self): return self._basis_card
    def knots (self): return self._vect     
    ###
    def is_open(self):
        #A knot vector is said to be open if 
        #its first and last knots appear p + 1 times.                       
        p = self.p()
        x = list(self._vect)
        if x.count(x[0]) != p+1 or x.count(x[-1]) != p+1:
            return False
        return True
    
    ###
    def xmin(self): return min(self._vect)
    def xmax(self): return max(self._vect)
                       
    ###
    def __len__(self): return len(self.knots())
    ###
    def show  (self): 
        print("polinomial degree : " , self.p())
        print("base caridnality  : " , self.n())
        print("knots             : " , self.knots())
    
def uniform_open_kv(xmin,xmax,p,n):
    v = np.concatenate((np.linspace(xmin,xmin,p),np.linspace(xmin,xmax,n-p+1),np.linspace(xmax,xmax,p)))
    return knot_vector(p,n,v)                     


# ## Bspline class

# In[17]:


import copy
import pandas as pd
#from scipy import integrate
#import itertools 
import time
#import pickle as pk
#import jsonpickle
#import json
import pickle
import re

class Bspline :
    
    ### constructor
    def __init__(self,
                 sh    = shape(),
                 knots = np.zeros(shape=(0,1),dtype=object),
                 cp    = None ,
                 prt   =  False
                ):
        
        # making sure knots is a np.ndarray
        #knots = np.asarray([knots])
        # call init_knot
        self  = self.init_knot(sh,knots)
        
        # decido se stampare a schermo tutto 
        self._print = prt
           
        #check control points dimension
        if len(self._cp.shape) != sh.dim() :
            print ("Control points error : wrong dimensions allocated")
            raise Exception()
            
        # assign control points
        if cp is not None :
            if cp.shape == self._cp.shape :
                self._cp = cp
            else :
                print ("Control points error : wrong dimensions")
                raise Exception()        
    ### function called by __init__
    def init_knot(self,sh,knots):
        
        # assign value to class members
        #dim    = sh.dim()        
        self._sh = sh
        
        #check dimension
        if len(knots) > sh.dim() :
            print("Warning : knots lenght too long")
        elif len(knots) < sh.dim() :
            print("Error : knots lenght too short")
            raise Exception()            
         
        #copy knot vectors
        self._kv = knots[0:sh.dim()]
        
        #print("len(self._kv) :",len(self._kv))
        #print("len(knots) :",len(knots))
            
        #allocate init 
        init     = self.Index_t()        
        
        # range over all "dimensions"
        for i in range(0,sh.dim()):
             
            # getting some values...
            p = knots[i].p()
            n = knots[i].n()
            k = knots[i].knots()
            v = len(k)

            # the i^th dimension will be n-long
            init[i] = n
            
            # checking everything is ok
            if not np.all(np.diff(k) >= 0) :
                raise Exception("knot vector not sorted")
            # va bene anche con p == 0
            if not p >= 0 : 
                raise Exception("wrong polyniomial degree : ",p)
            if not n > 1 : 
                raise Exception("wrong basis cardinality :",n)
            if not p+n+1 == v :
                raise Exception("wrong knot vector size :",p+n+1,"|=",v)                

        # assign control points dimension
        #print("init :" , init)
        # assign control points dimension
        #print("init :" , init)
        #self.clear_cp()
        self._cp = np.zeros(shape=tuple(map(int,init)),dtype=object)
        self._cp.fill(self.Type_out())
        if self.codim() == 1 :
            self._cp = self._cp.astype(float)
        
        #
        self._stiffness_matrix = None
        self._overlap_matrix = None
        #self._load_vector = None
        self._ready_sm = False #stiffness matrix
        self._ready_om = False #overlap matrix
        #self._ready_lv = False #load vector
        
        self._trace_Bspline = [ 0 for i in range(self.dim())]
        self._ready_trace = [ False for i in range(self.dim())]
        #print("self._cp : ",self._cp)
        return self    
    ###
    def copy(self):
        return copy.copy(self)    
    ### some function imitating C++ (template) type initialization
    def Type_out(self): #ok
        #if self.codim() == 1 :
        #    return 0.
        #else :
        #tenere così
        return np.zeros(self._sh.codim(),dtype=float)
    #def Type_in  (self) : #da rifare
    #    return np.zeros(shape=(self._sh.dim(),1))
    def Index_t  (self) : #ok
        return np.zeros(shape=(self._sh.dim(),1),dtype=int)    
    ### get control point value/coordinates
    def get_cp (self,index,check=True) :        
        try :
            # attenzione a passare degli int e non dei float
            return self._cp[index]
        except :
            if check == True :
                print("error : Bspline.get")
            return self.Type_out()         
    ### set control point value/coordinates
    def set_cp (self,index,value,check=True) :
        
        #if len(index) > 1 :
        #   print("error : Bspline.set, pass only one control point")
        #    raise Exception()
        #    return False
        
        #value = np.asarray(value)        
        #if len(value) == len(self._cp[index]) :
        value = np.asarray(value)
        value.reshape((self.codim(),1))
        if self.dim() != 1 :
            #print("convert to tuple:",type(index))
            ti = tuple(index)
            #print("type(ti):",type(ti))
            self._cp[ti] = value
            #print(self._cp)
        else :
            self._cp[index] = value
            #print("dim!=1")
            #print("index=",index)
            
        #print("index=",index)
        
        #else :
        #    if check == True :
        #        print("error : Bspline.set")
        #        raise Exception()
        #    return False
        self_ready_sm = False
        self_ready_om = False
        #self_ready_lv = False
        return True                                     
    ### some function returing some variables, not very usefull...
    def knots_vectors(self): return self._kv.copy()
    ###
    def get_knots_vector (self,index=None) :
        if index == None :        
            return self._kv
        else :
            return self._kv[index]
    ###
    def control_points (self, which="all") : 
        il = self.index_list(which)
        it = [ tuple(i) for i in il ]
        df = pd.DataFrame(index = it, columns=np.arange(0,self.codim()))
        for i in range(len(it)):
            index = df.index[i]
            df.iloc[i] = self._cp[index]
        return df
    ###
    def clear_cp (self) : 
        self._cp = np.zeros(self._cp.shape)
        self_ready_sm = False
        self_ready_om = False
        #self_ready_lv = False
    ###
    def show(self,what="all"):
        
        print("--------------------")
        
        #shape
        if what == "all" or "shape" in what :            
            self._sh.show()
            print("\n--------------------")
        
        #control points
        if what == "all" or "cp" in what :      
            print("control points")
            print("shape : ",self._cp.shape)
            print("coordinates")        
            print(self._cp)
            print("\n--------------------")
        
        #knot vectors
        if what == "all" or "kv" in what :    
            print("knot vectors")
            i=0
            for k in self._kv :
                print(i)
                print(k.show())
                i=i+1
            print("\n--------------------")        
    ###
    def dim(self)  : return self._sh.dim()
    ###
    def codim(self): return self._sh.codim()            
    ###
    def _find_k(self,x,t) :
        #print("x:",x)
        output = -1
        for i in range(0,len(t)-1) :
            if t[i] == t[i+1]:
                continue
            if t[i] <= x and x <= t[i+1] :
                output = i
        #if x >= t[-1] : #ultimo elemento
        #    output = -2
        return output
    ###
    def _deBoor(self,x,der=False):
        
        #
        if len(self._kv) > 1 :
            if self._print == True :
                print("deBoor error : knot vector size > 1")
            raise Exception()
            
        #
        t = self.get_knots_vector(0).knots()
        p = self.get_knots_vector(0).p()
        c = self._cp#self.control_points()
        
        #
        k = self._find_k(x,t)                
        #print("k:",k)
        if k < 0 :
            if self._print == True :
                print("deBoor error : k<0")
            #if k == -1 :
            return self.Type_out()
            #elif k == -2 :
            #    k = len(t)-2
            #    x = t[-1]
        
        #valuto la funzione
        #if der == False :
        return self._deBoor_private(k,x,t,c,p)
        #valuto la derivata prima della funzione
        #if der == True :
        #    return self._deBoor_private_derivative(k,x,t,c,p)        
    
    ### https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
    def _deBoor_private(self,k: int, x, t, c, p: int) :
        #Evaluates S(x).
        #
        #Arguments
        #---------
        #k: Index of knot interval that contains x.
        #x: Position.
        #t: Array of knot positions, needs to be padded as described above.
        #c: Array of control points.
        #p: Degree of B-spline.

        d = np.zeros(p+1,dtype=object)        
        for  j in range(0,p+1) :
            d[j] = c[j+k-p] if 0 <= j+k-p < len(c) else self.Type_out()
            
        
        for r in range(1, p + 1):            
            for j in range(p, r - 1, -1): 
                right = j+1+k-r
                left  = j+k-p
                if 0 <= left < len(t) and 0 <= right < len(t) :                    
                    alpha = (x - t[left]) / (t[right] - t[left])                    
                    d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]

        return d[p]
    ### https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/surface/bspline-de-boor.html
    def _iterative_deBoor(self,x) :
        
        # coordinate del punto in cui valutare la funzione
        # prendo tutte le coordinate tranne il primo
        x_surf = x[1:]      
        
        # knot vector
        # li prendo tutti tranne quello della prima dimensione
        curve_knot = self.knots_vectors()[1:]
        sub_kv     = self.get_knots_vector(0)
        new_sh     = shape(dim=self.dim()-1,codim=self.codim())        
        curve      = Bspline(new_sh,curve_knot)
        
        #
        m  = self.get_knots_vector(0).n()     
        #print("m:",m)
        #print("tuple:",tuple(map(int,[m])))
        surface_cp = np.zeros(shape=tuple(map(int,[m])),dtype=object)
        #surface_cp = np.zeros(shape=(m,1),dtype=object)
        surface_cp.fill(self.Type_out())
        #print("self.Type_out().shape:",self.Type_out().shape)
        #print("surface_cp[0].shape:",surface_cp[0].shape)
        
        #
        for i in range(0,m) :
            #print("curve._cp.shape:",curve._cp.shape)
            #print("self.get_cp(i).shape:",self.get_cp(i).shape)
            #print("\n")
            curve._cp = self.get_cp(i)      
            a = curve.evaluate(x_surf)
            #print("a:",a)
            #print("a.shape",a.shape)
            #print("len(surface_cp[i])",len(surface_cp[i]))
            surface_cp[i]=curve.evaluate(x_surf)
        
        #print("surface_cp : ",surface_cp)            
        out_sh = shape(dim=1,codim=self.codim())        
        return Bspline(out_sh,[sub_kv],surface_cp)    
    ###
    def evaluate(self,x):
        #sistemo dimensioni di x
        #if self._codim() > 1 :        
        X = self._correct_type_and_shape(x)
        #else :
        #    x = [X]
            
        #print("evaluate at :",X)
        if len(X) == 1 :
            x = X[0]
            #print("one x value passed : ",X," -> ",x)            
            # ho passato solo un valore
            if self.dim() == 1 :                
                out = self._deBoor(x)
                #print("dim == 1, out : ", out) 
                #return self._deBoor(x)
            else :
                #print("dim : ",self.dim())    
                curve_final = self._iterative_deBoor(x)
                #print("curve_final : ")
                #curve_final.show()
                out = curve_final._deBoor(x[0])
                #out = curve_final._deBoor(x[0])
                #print("dim >= 1, out : ", out) 
                #return curve_final._deBoor(x[0])
            if self.codim() == 1 :
                out = float(out)
            return out
        else :
            #print("x vector of lenght ",len(X)," passed")
            #if self.codim() == 1 :
            #    out = [float(self.evaluate(j)) for j in X ]
            #else :
            out = [ self.evaluate(j) for j in X ]
            #if self.codim() > 1 :
            out = np.asarray(out)
            if self.codim() == 1 :
                out = out.reshape((len(out,)))
            return out
    ###
    def derivative(self,axis=-1):        
        
        #http://public.vrac.iastate.edu/~oliver/courses/me625/week5b.pdf        
        #der.clear_cp()

        out = list()
        #derK = der.copy()
        
        # axis
        if axis == -1 :
            axis = [ i for i in range(self.dim()) ]
        if hasattr(axis, '__len__') == False :
            axis = [axis]
        
        for K in axis : #range(0,self.dim()):

            der_kv = list()
            der_sh = self._sh     
            for k in range(0,self.dim()):

                kv     = self._kv[k]
                kt     = kv.knots()

                if k != K :
                    der_kv.append(knot_vector(kv.p(),kv.n(),kt)) 
                else :     
                    der_kv.append(knot_vector(kv.p()-1,kv.n()-1,kt[1:-1]))         
            #der_kv = knot_vector(kv.p()-1,kv.n()-1,kt[1:-1]) # for i in self._kv]                
            derK = Bspline(der_sh,der_kv)

            #
            KV     = self._kv[K]
            KT     = KV.knots()
            P      = KV.p()

            #print("KT :",KT)

            #
            #derK = der.copy()
            N = list()
            Ntot = 1
            for k in range(0,self.dim()):
                if k != K :
                #prendo solo quelli lungo un asse
                    kv     = derK._kv[k]
                    kt     = kv.knots()
                    p      = kv.p()
                    #devo ciclare su tutte le altre dimensioni

                    N.append(np.arange(0,kv.n()))
                    Ntot = Ntot * kv.n()
                    #X.append(kt)

            w = list(np.meshgrid(*N))

            #print(Ntot)
            index = np.zeros(shape=(Ntot,self.dim()-1),dtype=int)
            #print(index)
            for i in range(0,len(w)):
                a = w[i].reshape((Ntot,))
                index[:,i] = a
            #print(index)

            for i in index:
                left  = [ i[kk] for kk in range(0,len(i)) if kk < K  ] #i[0:K]
                right = [ i[kk] for kk in range(0,len(i)) if kk >= K  ] #i[K:]
                #print("K :",K,"->",left,"-",right)

                #kv     = self._kv[K]
                #kt     = kv.knots()
                #p      = kv.p()

                for k in range(0,KV.n()-1):
                    ii  = list(left) + [k] + list(right)
                    iip = list(left) + [k+1] + list(right)

                    #ii  = ii.reverse()
                    #iip = iip.reverse()
                    #print(k,"-",ii)
                    #a = self._cp[[0,0,0]]
                    #b = self._cp[tuple(ii)]
                    cp = P*(self._cp[tuple(iip)] - self._cp[tuple(ii)]) / ( KT[k+P+1] - KT[k+1] ) 
                    #print("K :",K,"-> k:",k,"->",ii,"-> cp:",cp, "from ",self._cp[tuple(iip)],\
                    #      " and " ,self._cp[tuple(ii)] )
                    ii = ii[0] if self.dim() == 1 else ii
                    derK.set_cp(ii,cp) 
            out.append(derK)

        if self.dim() == 1 or len(axis) == 1 :
            return out[0]
        else :
            return out
    ###
    def _transpose(self,left,right):
        if left==right :
            return self
        
        #swap knot vector
        new_kv = self._kv
        new_kv[right],new_kv[left]  = new_kv[left],new_kv[right] 
        
        #control points
        new_cp = self._cp
        
        #"trasponiamo" un array multidimensionale
        sh = new_cp.shape
        new_sh = sh
        new_sh[right],new_sh[left] = new_sh[left],new_sh[right]        
        new_cp = new_cp.transpose(new_sh)
        
        #
        return Bspline(new_sh,new_kv,new_cp)
    ### 
    def _correct_type_and_shape(self,x):
        X = np.asarray(x)
        X = X.reshape((int(X.size/self.dim()),self.dim()))
        return X
    #Galerkin method
    ###
    def index_list(self,which="all"):
        N = list()
        Ntot = 1
        for k in range(0,self.dim()):
            kv     = self._kv[k]
            #kt     = kv.knots()
            #p      = kv.p()
            #devo ciclare su tutte le altre dimensioni

            N.append(np.arange(0,kv.n()))
            Ntot = Ntot * kv.n()
            #X.append(kt)

        w = list(np.meshgrid(*N))

        #print(Ntot)
        index = np.zeros(shape=(Ntot,self.dim()),dtype=int)
        #print(index)
        for i in range(0,len(w)):
            a = w[i].reshape((Ntot,))
            index[:,i] = a
        #print(index)
        if which == "all" :
            return index
        else:
            edge = self.edge()
            if which == "edge" :
                return index[ edge["edge"] == True ]
            elif which == "int" :
                return index[ edge["edge"] == False ]            
    
    ###
    def basis_range(self):
            
        basis_range = list()#np.zeros(self.dim(),dtype=object)
        for k in range(0,self.dim()):
            kv = self._kv[k]
            kt = kv.knots()
            n  = kv.n()
            p  = kv.p() 
            data = {"index":range(0,n),                   "i-min":range(0,n),                   "i-max":np.arange(0,n)+p+1,                   "x-min":kt[0:n],                   "x-max":kt[p+1:p+1+n]}
            br = pd.DataFrame(data)#np.array(shape=(n,3))
            basis_range.append(br)

        #creo gli indici
        index_index = [ ("index", i) for i in range(0,self.dim()) ]
        index_min   = [ ("min"  , i) for i in range(0,self.dim()) ]
        index_max   = [ ("max"  , i) for i in range(0,self.dim()) ]
        index = index_index + index_min + index_max

        #
        mi = pd.MultiIndex.from_tuples(index)
        df = pd.DataFrame(columns=index)
        df = df.reindex(columns=mi)

        #
        il = self.index_list() 
        it = [ tuple(j) for j in il ]

        #
        for k in range(0,self.dim()):
            df[("index",k)] = il[:,k]           
            appo =  basis_range[k]
            #print(appo)
            for j in range(0,len(appo)):
                #print(j)
                df.loc[ df[("index",k)] == appo.loc[j,"index"] , ("min",k) ] = appo.loc[j,"x-min"]
                df.loc[ df[("index",k)] == appo.loc[j,"index"] , ("max",k) ] = appo.loc[j,"x-max"]
            
        #if self.dim() == 1 :
        #    return df #basis_range[0],
        #else :
        
        df["ii"] = it
        df.set_index("ii",inplace=True)
        
        return df #basis_range,df
    
    ###
    def adjacency_matrix(self):
        
        #
        i = self.index_list()#df["index"]
        it = [ tuple(j) for j in i ]
        am = pd.DataFrame(index=it,columns=it)

        #
        br = self.basis_range()
        #br["ii"] = it
        #br.set_index("ii",inplace=True)
        #am

        n = am.shape[1]
        for i in range(0,n):
            r = am.index[i]
            c = am.columns[i]

            a = list()
            for k in range(0,self.dim()):
                a.append( [ br.at[r,("min",k)] ,br.at[r,("max",k)] ] )    

            am.at[r,c] = True
            for j in range(i+1,n):   
                #print(i,"-",j)
                c = am.columns[j]
                intersection = True
                for k in range(0,self.dim()):
                    b = [ br.at[c,("min",k)] ,br.at[c,("max",k)] ]  
                    o = overlaps( a[k] ,b )
                    if o is None :
                        intersection = False
                        break
                am.at[r,c] = intersection
                am.at[c,r] = intersection
                
        return am
    
    ###
    def basis_max_min(self,r,br=None):
        
        if br is None :
            br = self.basis_range()
        
        a = list()
        for k in range(0,self.dim()):
            a.append( [ br.at[r,("min",k)] , br.at[r,("max",k)] ] )    
            
        return a
    
    ###
    def basis_overlap(self,r,c,br=None):
        
        if br is None :
            br = self.basis_range()
        
        rmm = self.basis_max_min(r,br)
        cmm = self.basis_max_min(c,br)
        o = list()
        for k in range(0,self.dim()):
            o.append( overlaps( rmm[k] ,cmm[k] ) )
            
        return o
    
    ###
    def overlap_matrix(self,opts=None):
        
        opts = self.prepare_opts(opts)
        
        if opts["ready_om"] == False:
            self._ready_om = False
            
        if self._ready_om == True :
            return self._overlap_matrix
        
        #definisco la funzione da integrare
        #definisco la norma di Lebesgue

        #if opts["norm"] == "L1" :
        integrate = lambda *xx : left.evaluate(xx)*right.evaluate(xx)
        #elif opts["norm"] == "L2" :
        #    integrate = lambda *xx : np.power(left.evaluate(xx)*right.evaluate(xx),2.0)

        #mi serve solo per avere una classe costrutira correttamente
        #der = self.derivative()        
        br = self.basis_range()

        smd = list()#stifness matrix for derivatives

        X = np.zeros(self.dim())
        conta = 1
        for k in range(0,self.dim()):
            conta = conta * opts["delta"][k]
        Xintegration = np.zeros(shape=(conta,self.dim()))



        if opts["print"] == True : print("preparation",end="\r")
        #for k in range(0,self.dim()):

        #if opts["print"] == True : print("\ndimension :",k)
        #d = der[k]
        #d.clear_cp()
        scalar = self._scalar()
        left   = scalar.copy()
        right  = scalar.copy()
        left.clear_cp()
        right.clear_cp()

        #adjacency matrix
        am = self.adjacency_matrix()

        smd1D = am.copy()    
        smd1D[ smd1D == False ] = 0.0 #np.nan
        smd1D[ smd1D == True ] = 0.0 #np.nan

        #calcolo il prodotto scalare dei gradienti
        n = am.shape[0]
        for i in range(0,n):
            r = am.index[i] 

            #creo la funzione di base
            #left.clear_cp()
            left.set_cp(r,1.0)

            for j in range(i,n):

                c = am.columns[j]

                if am.at[r,c] is False :
                    continue

                #print(i,"-",j)

                #creo la funzione di base
                #right.clear_cp()
                right.set_cp(c,1.0)
                ov = self.basis_overlap(r,c,br)

                #
                X = [ np.delete(np.linspace(ov[k][0],ov[k][1],                                            opts["delta"][k]+1,endpoint=False),0)                      for k in range(0,self.dim()) ]

                area = 1
                for k in range(0,self.dim()):
                    area = area * ( ov[k][1] - ov[k][0] )

                m = np.meshgrid(*X)
                for k in range(0,self.dim()):
                    Xintegration[:,k] = m[k].flatten()

                #
                if opts["print"] == True : start = time.time()

                y = integrate (Xintegration)
                #print(y)
                #uso sempre L1
                #if opts["norm"] == "L1" :
                res = np.mean(y)*area
                #elif opts["norm"] == "L2" :
                #print("L2 norm")
                #res = None #np.sqrt(np.mean(y))*area

                #res = integrate.nquad( integral , ov , opts = opts)[0] #solo risultato
                if opts["print"] == True : endt = time.time()
                if opts["print"] == True : print(r,"-",c," -> ", endt - start," s")

                #    
                smd1D.at[r,c] = res
                smd1D.at[c,r] = res #matrice simmetrica

                #cancello
                right.set_cp(c,0.0)

                #cancello
            left.set_cp(r,0.0)

        self._overlap_matrix = smd1D.copy()
        self._ready_om = True
        return smd1D
            #smd.append(smd1D)
        
    ###
    def prepare_opts(self,opts,opts2=None):
        if opts is None:
            opts = {}
        if "print" not in opts :
            opts["print"] = False
        if "delta" not in opts :
            opts["delta"] = np.full((self.dim(),),4)
        if "norm" not in opts :
            opts["norm"] = "L1"
        if "del-edge" not in opts :
            opts["del-edge"] = True
        if opts2 is not None :
            opts.update(opts2)
        if "ready_om" not in opts :
            opts["ready_om"] = True
        if "ready_sm" not in opts :
            opts["ready_sm"] = True
        #if "ready_lv" not in opts :
        #    opts["ready_lv"] = False
        if "ready_trace" not in opts :
            opts["ready_trace"] = [ False for i in range(self.dim())]
        
        return opts
    
    ###
    def edge(self):
        il = self.index_list()
        df = pd.DataFrame( il , index = tuple(il) ,columns = np.arange(0,self.dim()) )
        df["edge"] = True
        #df["corner"] = True
        allx = np.arange(0,self.dim())
        e = np.zeros(self.dim())
        for i in range(0,df.shape[0]):
            for k in allx :                
                kv = self._kv[k]
                e[k] = df.iloc[i,k] <= 0 or df.iloc[i,k] >= kv.n()-1
                #e[k] = df.iloc[i,k] < kv.p() or df.iloc[i,k] > kv.n()-kv.p()-1
            df.at[df.index[i],"edge"] = np.any(e)
            #df.at[df.index[i],"corner"] = np.all(e)
        df.drop(columns=allx,inplace=True)
        return df       
    
    ###
    def approximate(self,func,opts=None):
        #http://hplgit.github.io/INF5620/doc/pub/sphinx-fem/._main_fem002.html
        opts = self.prepare_opts(opts)

        om = self.overlap_matrix(opts)
        lv = self.load_vector(func,opts)

        om.replace(np.nan,0.0,inplace=True)
        omnp = np.asarray(om)

        if self.codim() == 1 :

            lv.replace(np.nan,0.0,inplace=True)

            lvnp = np.asarray(lv)

            cp  = np.linalg.solve(omnp,lvnp)
            out = pd.DataFrame(cp,index=om.index,columns=["cp"])

            for i in range(len(cp)):
                j = out.index[i]
                self._cp[j]  = out.at[j,"cp"]

        else :

            lv2 = pd.DataFrame(columns=np.arange(0,self.codim()),index=lv.index)
            for k in range(self.codim()):
                for i in lv2.index :
                    lv2.at[i,k] =lv.at[i,"cp"][k]
            lv2.replace(np.nan,0.0,inplace=True)

            out = pd.DataFrame(index=om.index,columns=np.arange(0,self.codim()))

            for k in range(self.codim()):

                lvnp = np.asarray(lv2[k])
                cp  = np.linalg.solve(omnp,lvnp)

                out[k] = cp

            for i in range(len(cp)):
                j = out.index[i]
                self.set_cp(j, out.iloc[i])
        if self.codim() == 1 :
            self._cp = self._cp.astype(float)
        return out  
       
    ###
    def _scalar(self):#restituisce la stessa Bspline ma con cp scalari
        sh = shape(dim=self.dim(),codim=1)
        kv = self._kv        
        return Bspline(sh,kv)#cp: default value
    
    ###
    def get_gD_shape(self):
        index = self.index_list("edge")
        it = [ tuple(j) for j in index ]
        columns = ["Dirichlet"]
        return pd.DataFrame(data=np.full(len(it),0),index=it,columns=columns)
        
    ###
    def _get_xint(self):
        
        # restituisce i punti di bordo interpolatori
        
        #index list
        index_list = self.index_list("edge")
        il = pd.DataFrame(index_list)
        #index list (tuple)
        it =  [ tuple(i) for i in index_list ]

        #xmin xmax
        xmin_vect = [ self._kv[i].xmin() for i in range(self.dim()) ]
        xmax_vect = [ self._kv[i].xmax() for i in range(self.dim()) ]

        # x interpolatorio
        xint = pd.DataFrame(np.zeros(shape=(len(il),self.dim())),index = it)

        #percorro la matrice "il" dall'alto verso il basso, una colonna alla volta
        for k in range(self.dim()): #ciclo sulle dimensioni/componenti dei punti : COLONNE
            n = self._kv[k].n()
            xmin = xmin_vect[k]
            xmax = xmax_vect[k]
            knot0 = self.get_knots_vector(k).knots()
            p = self.get_knots_vector(k).p()
            n = self.get_knots_vector(k).n()
            knot  = knot0[p:n+1]
            for i in range(len(il)) : #ciclo sui punti di bordo : RIGHE
                if il.at[i,k] == 0 :
                    xint.iloc[i,k] = xmin
                elif il.at[i,k] == n-1 :
                    xint.iloc[i,k] = xmax
                else :
                    xint.iloc[i,k] = knot[xint.index[i][k]]#knot[k]
                    
        return xint
    
    ###
    def trace(self,n,opts=None):
        
        
        opts = self.prepare_opts(opts)
        
        if opts["ready_trace"][n] == False:
            self._ready_trace[n] = False
            
        if self._ready_trace[n] == True :
            return self._trace_Bspline[n]
        
        # operazione di traccia: 
        # https://en.wikipedia.org/wiki/Trace_operator
        # restringo il dominio di una funzione al suo bordo
        # una superficie diventa una curva
        
        # n : dimensione da tracciare
        
        new_sh = shape(dim=self.dim()-1,codim=self.codim())        
        kv     =    [ self.knots_vectors()[ii] for ii in range(0,self.dim()) if ii != n ]    
        curve  = Bspline(new_sh,kv)
        
        self._ready_trace[n] = True
        self._trace_Bspline[n] = curve.copy()
        
        return curve
        
    
    ### 
    def Dirichlet_BC(self,gD,opts):
        
        # to_approx: 
        def to_approx_private(xx,x_min_max,kk,gD):        
            xyz = np.zeros(shape=(len(xx),self.dim()))
            jj = 0
            for ii in range(self.dim()):
                if ii != kk :
                    xyz[:,ii] = xx[:,jj]
                    jj = jj+1
                else :
                    xyz[:,kk] = np.full(len(xx),x_min_max) 
            return gD(xyz)
        
        opts = self.prepare_opts(opts)

        new_sh     = shape(dim=self.dim()-1,codim=self.codim())        
        for k in range(self.dim()):
            #k = i 
            #kv = [ self.knots_vectors()[ii] for ii in range(0,self.dim()) if ii != k ]    
            #curve      = Bspline(new_sh,kv)
            curve = self.trace(k,opts)
            #
            i1 = 0
            i2 = self.knots_vectors()[k].n()-1
            xmin = min(self.knots_vectors()[k].knots())
            xmax = max(self.knots_vectors()[k].knots())
            # xmin,xmax
            for index,x_min_max in zip([i1,i2],[xmin,xmax]) :

                xmin = min(self.knots_vectors()[k].knots())
                to_approx = lambda xx : to_approx_private(xx,x_min_max,k,gD)
                curve.approximate(to_approx,opts)

                cp = curve.control_points()
                for i in cp.index :
                    j = list(i)
                    j.insert(k,index)
                    self.set_cp(j,cp[cp.index == i])

        cp = self.control_points()
        ii = self.index_list("edge")
        it = [ tuple(i) for i in ii ]
        return cp.loc[it]

        #xint = self._get_xint()        
        # convert into numpy array
        #xintnp = np.array(xint)        
        # calcolo gD nei punti interpolatori
        #yintnp = gD(xint)        
        #converto in DataFrame
        #yint = pd.DataFrame(yintnp,index=xint.index)        
        #return yint
        
    ###
    def Galerkin(self,f,gD=None,opts=None):
        
        #               u = unknown function
        #          -Lap u = f    on Omega
        #               u = gD   on Dirichlet boundary
        
        opts = self.prepare_opts(opts)

        #edge
        edge = self.edge()

        #stiffness matrix
        sm = self.stiffness_matrix(opts) 
        sm.replace(np.nan,0.0,inplace=True)

        #indici dei dof interni e di bordo
        index_int  = edge.index[ edge["edge"] == False ]
        index_edge = edge.index[ edge["edge"] == True  ]

        #degrees of freedem: internal
        dof_int = sm.copy()
        dof_int.drop( index_edge ,inplace = True,axis=0)
        dof_int.drop( index_edge ,inplace = True,axis=1)

        #degrees of freedem: edge
        dof_edge = sm.copy()
        dof_edge.drop( index_edge ,inplace = True,axis=0)
        dof_edge.drop( index_int  ,inplace = True,axis=1)

        #load vector
        lv = self.load_vector(f,opts)        
        lv.replace(np.nan,0.0,inplace=True)
        lv.drop( index_edge  ,inplace = True)  

        #convert into numpy array
        dinp = np.asarray(dof_int)
        denp = np.asarray(dof_edge)
        #
        if gD is None :
            gD = lambda xx : np.zeros(shape=(len(xx),self.codim()))
        
        # valore di gD nei punti interpolatori
        # ERRORE:
        # la Bspline non è iterpolatorio solo negli angoli
        # devo creare una Bspline di dimensione dim-1
        # approssimare la funzione gD tramite questa Bspline
        # i control points di questa diventeranno gDv
        # cioè i control points di bordo della Bspline originaria
        gDv = self.Dirichlet_BC(gD,opts)
        
        #load vector
        lvnp = np.asarray(lv)
        #punti di bordo
        gDnp = np.asarray(gDv,dtype=float)
        #prodotto righe per colonne
        edlv = np.dot(denp,gDnp)
        
        #solve linear system
        cpint = np.linalg.solve(dinp,lvnp-edlv) 

        #preparo gli indice della variabile di output
        index = self.index_list()
        it = [ tuple(j) for j in index ]
        out = pd.DataFrame(index=it,columns=["cp"])

        #assegno ai control points i valori calcolati
        #valori interni
        for i in range(len(index_int)):
            j = index_int[i]
            self._cp[j] = cpint[i]
            out.iloc[out.index == j] = cpint[i]
            
        #assegno ai control points i valori calcolati
        #valori di bordo interpolanti
        for i in range(len(index_edge)):
            j = index_edge[i]
            self._cp[j] = gDnp[i]
            out.iloc[out.index == j] = gDnp[i]
        #
        if self.codim() == 1 :
            self._cp = self._cp.astype(float)
        return out      
    
    ###
    def stiffness_matrix(self,opts=None):
        
        opts = self.prepare_opts(opts)
        
        # controllo se ho già calcolato la matrici di stiffness
        if opts["ready_sm"] == False:
            self._ready_sm = False
            
        if self._ready_sm == True :
            return self._stiffness_matrix

        # omd : overlap matrix of derivatives
        omd = list()
        der = self.derivative() if self.dim() > 1 else [self.derivative()]
        if opts["print"] == True : print("preparation",end="\r")
        for k in range(0,self.dim()):
            if opts["print"] == True : print("\ndimension :",k)  
            omd1D = der[k].overlap_matrix(opts)
            omd.append(omd1D)
        del omd1D
        # smd : matrice di overlap delle derivate

        ###
        # am: adjacency matrix
        am = self.adjacency_matrix()
        # sm1D : stiffness matrix 1D
        # considero le derivate parziali lungo solo un asse
        sm1D = pd.DataFrame(0.0,index=am.index,columns=am.columns)
        n = sm1D.shape[0]#quadrata

        # escludo subito dal calcolo i termini di bordo 
        # impostando a nan il valore in sm1D
        # edge : questo non dovrebbe servire più
        #edge = self.edge()


        #sm1D = am.copy()    
        #sm1D[ sm1D == False ] = 0.0 #None

        #stiffness matrix: sum over all dimensione
        #matrice stiffness finale, è la somma della matrici parziali (su una sola dimensione)
        #questa è in realtà una lista contenente le matrici parziali

        # sm : stiffness matrix
        # lista che contiene tutte le sm1D
        sm = list()



        # creo due copie della Bspline
        left  = self.copy()
        left.clear_cp()
        right = self.copy()
        right.clear_cp()


        if opts["print"] == True : print("\nstifness matrix",end="\r")
        for k in range(0,self.dim()):
            if opts["print"] == True : print("\ndimension :",k)

            #ciclo su tutte le funzioni di base
            for i in range(0,n):

                #indice della funzione di base di sinista (r=righe)
                r = am.index[i] 

                #controllo che il termine non sia di bordo
                #if edge.at[r,"edge"] == True :
                #    sm1D.at[r,:] = np.nan
                #    sm1D.at[:,r] = np.nan
                #    continue

                #left.clear_cp()
                left.set_cp(r,1.0)
                # calcolo la derivata della funzione di base di sinistra
                # dl : derivatives left
                # ho modificato gli input della funzione derivative
                #dl = left.derivative()[k] if self.dim() > 1 else left.derivative()
                dl = left.derivative(axis=k) #if self.dim() > 1 else left.derivative()

                #cancello, tanto non mi serve valutarla
                #mi servono soltanto i coefficienti della derivata
                left.set_cp(r,0.0)


                # ATTENZIONE :
                # qui inizia l'algoritmo vero e proprio


                #control points, modifico il tipo e cerco quelli non nulli
                #cpl = dl._cp.astype(float) 
                # nzcpl : non zero control points (left) indices
                #nzcpl = np.argwhere( cpl != 0.).tolist() 

                cpl = dl.control_points()
                cpl.columns = ["value"]
                nzcpl = cpl[ cpl["value"] !=0. ]


                #ciclo su tutte le altre funzioni di base
                for j in range(i,n):
                    #if opts["print"] == True : print(i,"-",j,end="\r")

                    #indice della funzione di base di destra (c=colonne)
                    c = am.columns[j]   

                    #controllo che il termine non sia di bordo
                    #if edge.at[c,"edge"] == True :
                    #    sm1D.at[c,:] = np.nan
                    #    sm1D.at[:,c] = np.nan
                    #    continue


                    #right.clear_cp()
                    right.set_cp(c,1.0)
                    # ho modificato gli input della funzione derivative
                    #dr = right.derivative()[k] if self.dim() > 1 else right.derivative()
                    dr = right.derivative(k) #if self.dim() > 1 else right.derivative()
                    #cancello, tanto non mi serve valutarla
                    #mi servono soltanto i coefficienti della derivata
                    right.set_cp(c,0.0)

                    #control points, modifico il tipo e cerco quelli non nulli
                    #cpr = dr._cp.astype(float)
                    # nzcpr : non zero control points (right) indices
                    #nzcpr = np.argwhere( cpr != 0.).tolist() 

                    cpr = dr.control_points()
                    cpr.columns = ["value"]
                    nzcpr = cpr[ cpr["value"] !=0. ]

                    #attenzione all'ordine
                    #genero tutte le coppie
                    #allpairs = list(itertools.product(nzcpl,nzcpr))
                    if opts["print"] == True : print(r,"-",c)
                    if opts["print"] == True : print(cpl,"->",nzcpl)
                    if opts["print"] == True : print(cpr,"->",nzcpr)

                    for rr in nzcpr.index:
                        #right value
                        rv = cpr.at[rr,"value"] 
                        for ll in nzcpl.index:
                            #left value
                            lv = cpl.at[ll,"value"]  
                            #derivative value
                            dv = omd[k].at[ll,rr]
                            #
                            sm1D.at[r,c] = sm1D.at[r,c] + lv*rv*dv

                    #sm1D.at[r,c] = 0.0         
                    #for w in range(0,len(allpairs)):
                    #    li = tuple(allpairs[w][0])
                    #    ri = tuple(allpairs[w][1])
                    #    ll = cpl[li]
                    #    rr = cpr[ri]                
                    #    #if ll != 0.0 and rr != 0.0 :                    
                    #    dd = omd[k].at[li,ri]
                    #   if opts["print"] == True : print("sum : ",ll," | ",rr,"|",dd)
                    #    #if dd is None : dd = 0.0
                    #    sm1D.at[r,c] = sm1D.at[r,c] + ll*rr*dd

                    sm1D.at[c,r] = sm1D.at[r,c]
                    if opts["print"] == True : print("\n")

            sm.append(sm1D)
            out = sum(sm)
            #
            #out.drop(edge.index[ edge["edge"] == True ],inplace = True,axis=0)
            #out.drop(edge.index[ edge["edge"] == True ],inplace = True,axis=1)

        self._stiffness_matrix = out.copy()
        self._ready_sm = True
        return out
     
    ###
    def load_vector(self,func,opts=None): 
        
        opts = self.prepare_opts(opts)
        
        # controllo se ho già calcolato la matrici di stiffness
        #if opts["ready_lv"] == False:
        #    self._ready_lv = False
            
        #if self._ready_lv == True :
        #    return self._load_vector

        # escludo subito dal calcolo i termini di bordo 
        #impostando a nan il valore in sm1D
        #if opts["del-edge"] == True :
        edge = self.edge()

        #
        X = np.zeros(self.dim())
        conta = 1
        for k in range(0,self.dim()):
            conta = conta * opts["delta"][k]
        Xintegration = np.zeros(shape=(conta,self.dim()))

        #
        br = self.basis_range()
        il = self.index_list()
        #if self.dim() > 1 :
        il = [ tuple(i) for i in il ]
        #else :
        #    il = [ int(i) for i in il ]
        scalar = self._scalar()
        scalar.clear_cp()
        #piccola modifica : self -> scalar
        out = pd.DataFrame(columns=["index","cp"],dtype=object)
        out["index"] = il
        #out = pd.DataFrame({"index" : il , "cp" : scalar._cp.copy().flatten().astype(float) })
        out.set_index("index",inplace=True)

        #definisco la norma di Lebesgue

        def integrate_ND(xx):
            a = scalar.evaluate(xx)
            b = func(xx)
            out = b.copy()
            for i in range(self.codim()):
                out[:,i] = a*b[:,i]
            return out

        def integrate_1D(xx):
            A = scalar.evaluate(xx)
            B = func(xx)
            return [ float(a*b) for a,b in zip(A,B) ] 

        if self.codim() == 1 :
            integrate = integrate_1D
        else :
            integrate = integrate_ND
            #print("codim > 1: Galerkin method is not defined")
            #raise Exception()
        #integrate = lambda *xx : 


        for i in il :

            #controllo che il termine non sia di bordo
            #if edge.at[i,"edge"] == True :
            #    out.at[i,"cp"] = np.nan
            #    continue

            #i = tuple(i)
            scalar.set_cp(i,1.0)
            ov = self.basis_overlap(i,i,br) #overlap
            
            # ho messo endpoint = True
            # rimetto endpoint = False
            # tenere assolutamente endpoint = False
            X = [ np.delete(np.linspace(ov[k][0],ov[k][1],opts["delta"][k]+1,endpoint=False),0)                  for k in range(0,self.dim()) ]
            area = 1
            for k in range(0,self.dim()):
                area = area * ( ov[k][1] - ov[k][0] )
            m = np.meshgrid(*X)
            for k in range(0,self.dim()):
                Xintegration[:,k] = m[k].flatten()

            #
            if opts["print"] == True : start = time.time()
            y = integrate (Xintegration)
            if opts["norm"] == "L1" :
                #modifica
                res = np.mean(y, axis=0) if self.codim() > 1 else np.mean(y)
            elif opts["norm"] == "L2" :
                #modifica
                res = np.sqrt(np.sum(np.power(y,2.0),axis=0))/len(y) if self.codim() > 1                 else np.sqrt(np.sum(np.power(y,2.0)))/len(y)
            #modifica
            out.at[i,"cp"] = res * area
            if opts["print"] == True : endt = time.time()
            if opts["print"] == True : print(i," -> ", endt - start," s")                
            #out[i] = integrate.nquad( func , ov , opts = opts)[0] 

            #cancello
            scalar.set_cp(i,0.0)
        #if opts["del-edge"] == True:
        #    out = out.drop(edge.index[ edge["edge"] == True ])  
            
        #self._load_vector = out.copy()
        #self._ready_lv = True
        return out       
    
    
    ###
    def save(self,variable,filename):
        if variable == "sm":
            var = self._stiffness_matrix
        elif variable == "om":
            var = self._overlap_matrix
        elif variable == "lv":
            var = self._load_vector
        elif variable == "cp" :
            var = self._cp
        var.to_csv(filename,index_label="index")
        
    ###
    def load(self,variable,filename):
        
        var = pd.read_csv(filename)
        var.index = [tuple_from_str(i) for i in var.iloc[:,0]]
        var = var.drop('index',axis=1)
        
        if variable == "sm":
            var.columns = [tuple_from_str(i) for i in var.columns]
            self._stiffness_matrix = var.copy()
            self._ready_sm = True
        elif variable == "om":
            var.columns = [tuple_from_str(i) for i in var.columns]
            self._overlap_matrix = var.copy()
            self._ready_om = True
        elif variable == "lv":
            self._load_vector = var.copy()
            self._ready_lv = True
        elif variable == "cp" :
            index = var.index
            for i in index:
                self._cp[i] = var.at[i,"cp"]
            #self._cp = var.copy()


###
def tuple_from_str(string):
    return tuple(map(int, re.findall(r'[0-9]+', string)))

    
###        
def overlaps(a, b):
    c = [ max(a[0],b[0]) , min(a[1],b[1]) ]
    
    if c[0] >= c[1] :
        return None
    else :
        return c


# 
#         #edge
#         edge = self.edge()
# 
#         #indici dei dof interni e di bordo
#         index_int  = edge.index[ edge["corner"] == False ]
#         index_edge = edge.index[ edge["corner"] == True  ]        
#         
#         #overlap matrix   
#         om = self.overlap_matrix(opts)
#         om.replace(np.nan,0.0,inplace=True)
# 
#         #degrees of freedem: internal
#         dof_int = om.copy()
#         dof_int.drop( index_edge ,inplace = True,axis=0)
#         dof_int.drop( index_edge ,inplace = True,axis=1)
# 
#         #degrees of freedem: edge            
#         dof_edge = om.copy()
#         dof_edge.drop( index_edge ,inplace = True,axis=0)
#         dof_edge.drop( index_int  ,inplace = True,axis=1)
# 
#         #convert into numpy array
#         dinp = np.asarray(dof_int)
#         denp = np.asarray(dof_edge)
#         
#         #load vector
#         lv = self.load_vector(func,opts)
#         lv.drop( index_edge  ,inplace = True) 
#         lvnp = np.asarray(lv["cp"]).astype(float)
# 
#         #
#         XY = pd.DataFrame(index = index_edge,columns=np.arange(0,self.dim()))
#         for k in range(self.dim()):
#             kv = self.knots_vectors()[k]
#             nmax = kv.n()-1
#             xmin = min(kv.knots())
#             xmax = max(kv.knots())
#             for i in index_edge:
#                 if i[k] == nmax :
#                     XY.at[i,k] = xmax
#                 else :
#                     XY.at[i,k] = xmin
# 
#         xy = np.asarray(XY).astype(float)
# 
#         gDnp = func(xy)
# 
#         #prodotto righe per colonne
#         edlv = np.dot(denp,gDnp)
# 
#         cpint = np.linalg.solve(dinp,lvnp-edlv) 
# 
#         #out = pd.DataFrame(cp,index=om.index,columns=["cp"])
# 
#         index = self.index_list()
#         it = [ tuple(j) for j in index ]
#         out = pd.DataFrame(index=it,columns=["cp"])
# 
#         for i in range(len(index_int)):
#             j = index_int[i]
#             self._cp[j] = cpint[i]
#             out.iloc[out.index == j] = cpint[i]
# 
#         #assegno ai control points i valori calcolati
#         #valori di bordo interpolanti
#         for i in range(len(index_edge)):
#             j = index_edge[i]
#             self._cp[j] = gDnp[i]
#             out.iloc[out.index == j] = gDnp[i]
#            
#         return out       
#        

# opts = self.prepare_opts(opts)        
# 
#         #edge
#         edge = self.edge()
# 
#         #indici dei dof interni e di bordo
#         index_int  = edge.index#[ edge["corner"] == False ]
#         index_edge = []#edge.index[ edge["corner"] == True  ]        
#         
#         #overlap matrix   
#         om = self.overlap_matrix(opts)
#         om.replace(np.nan,0.0,inplace=True)
# 
#         #degrees of freedem: internal
#         dof_int = om.copy()
#         dof_int.drop( index_edge ,inplace = True,axis=0)
#         dof_int.drop( index_edge ,inplace = True,axis=1)
# 
#         #degrees of freedem: edge            
#         dof_edge = om.copy()
#         dof_edge.drop( index_edge ,inplace = True,axis=0)
#         dof_edge.drop( index_int  ,inplace = True,axis=1)
# 
#         #convert into numpy array
#         dinp = np.asarray(dof_int)
#         denp = np.asarray(dof_edge)
#         
#         #load vector
#         lv = self.load_vector(func,opts)
#         lv.drop( index_edge  ,inplace = True) 
#         lvnp = np.asarray(lv["cp"]).astype(float)
# 
#         #
#         XY = pd.DataFrame(index = index_edge,columns=np.arange(0,self.dim()))
#         for k in range(self.dim()):
#             kv = self.knots_vectors()[k]
#             nmax = kv.n()-1
#             xmin = min(kv.knots())
#             xmax = max(kv.knots())
#             for i in index_edge:
#                 if i[k] == nmax :
#                     XY.at[i,k] = xmax
#                 else :
#                     XY.at[i,k] = xmin
# 
#         xy = np.asarray(XY).astype(float)
# 
#         gDnp = func(xy)
# 
#         #prodotto righe per colonne
#         edlv = np.dot(denp,gDnp)
# 
#         cpint = np.linalg.solve(dinp,lvnp-edlv) 
# 
#         #out = pd.DataFrame(cp,index=om.index,columns=["cp"])
# 
#         index = self.index_list()
#         it = [ tuple(j) for j in index ]
#         out = pd.DataFrame(index=it,columns=["cp"])
# 
#         for i in range(len(index_int)):
#             j = index_int[i]
#             self._cp[j] = cpint[i]
#             out.iloc[out.index == j] = cpint[i]
# 
#         #assegno ai control points i valori calcolati
#         #valori di bordo interpolanti
#         for i in range(len(index_edge)):
#             j = index_edge[i]
#             self._cp[j] = gDnp[i]
#             out.iloc[out.index == j] = gDnp[i]

# ### Derivative

# $\mathbf{r}\left(t\right) = \sum_{k=0}^{n} \mathbf{P}_k N^{p}_{k}\left(t\right)$
# 
# $\dfrac{ \partial \mathbf{r}}{\partial t} =\sum_{k=0}^{n} \mathbf{P}_k \frac{ \partial  N^{p}_{k}\left(t\right) }{\partial t}$
# 
# $ N^{p}_{k}\left(t\right) = \dfrac{t - u_k}{u_{k+p}-u_{k}} N^{p-1}_{k}\left(t\right) - \dfrac{t - u_{k+p+1}}{u_{k+p+1}-u_{k+1}} N^{p-1}_{k+1}\left(t\right) $
# 
# Se $p=1$ otteniamo
# $\dfrac{ \partial  N^{p}_{k}\left(t\right) }{\partial t} = p \left[ \dfrac{N^{p-1}_{k}}{u_{k+p}-u_{k}}  - \dfrac{N^{p-1}_{k+1}}{u_{k+p+1}-u_{k+1}}  \right] $
# 
# Procediamo per induzione (tenendo conto che passiamo da polinomi di grado $p$ a polinomi di grado $p-1$ e che stiamo togliendo il primo e ultimo elemento del knot vector)
# 
# $\dfrac{ \partial  N^{p}_{k} }{\partial t} =  %
# \dfrac{N^{p-1}_{k}}{u_{k+p}-u_k} + %
# \dfrac{t-u_k}{u_{k+p}-u_k}\dfrac{ \partial  N^{p-1}_{k} }{\partial t} - %
# \dfrac{N^{p-1}_{k+1}}{u_{k+p+1}-u_{k+1}} - %
# \dfrac{t-u_{k+1}}{u_{k+p+1}-u_{k+1}}\dfrac{ \partial  N^{p-1}_{k+1}}{\partial t}  \\%
# = \dfrac{N^{p-1}_{k}}{u_{k+p}-u_k} - %
# \dfrac{N^{p-1}_{k+1}}{u_{k+p+1}-u_{k+1}} + %
# \dfrac{\left(t-u_k\right)\left(p-1\right)}{u_{k+p}-u_k} \left[ \dfrac{N^{p-2}_{k}}{u_{k+p}-u_{k}}  - \dfrac{N^{p-2}_{k+1}}{u_{k+p+1}-u_{k+1}} \right] - %
# \dfrac{\left(t-u_{k+1}\right)\left(p-1\right)}{u_{k+p+1}-u_{k+1}} \left[ \dfrac{N^{p-2}_{k+1}}{u_{k+p+1}-u_{k+1}}  - \dfrac{N^{p-2}_{k+2}}{u_{k+p+2}-u_{k+2}} \right] \\
# = \dfrac{N^{p-1}_{k}}{u_{k+p}-u_k} - %
# \dfrac{N^{p-1}_{k+1}}{u_{k+p+1}-u_{k+1}} + %
# \dfrac{\left(p-1\right)}{u_{k+p}-u_k} N^{p-1}_{k} -%
# \dfrac{\left(p-1\right)}{u_{k+p+1}-u_{k+1}}N^{p-1}_{k+1} \\
# = p \left[ \dfrac{N^{p-1}_{k}}{u_{k+p}-u_k} - %
# \dfrac{N^{p-1}_{k+1}}{u_{k+p+1}-u_{k+1}} \right]$
# 
# dunque
# 
# $\dfrac{ \partial \mathbf{r}}{\partial t} = %
# \sum_{k=0}^{n} \mathbf{P}_k p \left[ \dfrac{N^{p-1}_{k}}{u_{k+p}-u_k} - %
# \dfrac{N^{p-1}_{k+1}}{u_{k+p+1}-u_{k+1}} \right]  \\
# = \sum_{k=0}^{n-1} N^{p-1}_{k+1} \left[ p \dfrac{\mathbf{P}_{k+1} - \mathbf{P}_{k}}{u_{k+p+1}-u_{k+1}}\right] \\
# = \sum_{k=0}^{n-1} N^{p-1}_{k} \mathbf{Q}_{k} $
# 
# poiché $N^{p-1}_{k+1}$ valutato sul knot vector originale è uguale a $N^{p-1}_{k}$ valutato sul nuovo (open) knot vector.
# Abbiamo definito
# 
# $\mathbf{Q}_{k} = p \dfrac{\mathbf{P}_{k+1} - \mathbf{P}_{k}}{u_{k+p+1}-u_{k+1}} $

# ### Surface Bspline

# $\mathbf{r}\left(x,y\right) = \sum_{i=0}^{n}\sum_{j=0}^{m} \mathbf{P}_{ij} N^{p}_{i}\left(x\right) M^{q}_{j}\left(y\right) \\
# = \sum_{i=0}^{n} N^{p}_{i}\left(x\right) \left( \sum_{j=0}^{m} \mathbf{P}_{ij} M^{q}_{j}\left(y\right) \right)  \\ 
# = \sum_{i=0}^{n} N^{p}_{i}\left(x\right) \mathbf{Q}_{i}\left(y\right)$
# 
# con 
# 
# $ \sum_{j=0}^{m} \mathbf{P}_{ij} M^{q}_{j}\left(y\right) = \mathbf{Q}_{i}\left(y\right) $

# ## Galerkin method

# 
# 
