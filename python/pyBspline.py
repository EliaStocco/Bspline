#!/usr/bin/env python
# coding: utf-8

# In[4]:


get_ipython().system('jupyter nbconvert --to script pyBspline.ipynb')


# ## Shape class

# In[2]:


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

# In[3]:


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

def periodic_kv(xmin,xmax,p,n):
    #https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-closed.html
    v0 = np.linspace(0.,1.0,n+2*p+2,endpoint=True)
    v = xmin + (v0 - v0[p])/v0[n+1]*(xmax-xmin)
    return knot_vector(p,n+1+p,v,check_open=False) 

    #v0 = np.linspace(0.,1.0,n+2*p+2,endpoint=True)
    #v = xmin + (v0 - v0[p])/v0[n+1]*(xmax-xmin)
    #return knot_vector(p,n+p+1,v,check_open=False) 
    #v0 = np.linspace(0.,1.0,n+p+p,endpoint=True)
    #v = xmin + (v0 - v0[p])/v0[n]*(xmax-xmin)
    #return knot_vector(p,n+p+1,v,check_open=False)  
    #v0 = np.linspace(0.,1.0,n+1+p,endpoint=True)
    #v = xmin + (v0 - v0[p])/v0[n-p]*(xmax-xmin)
    #return knot_vector(p,n,v,check_open=False)  


# ## Bspline class

# In[2]:


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
from numpy import random
import scipy.special

class Bspline :
    
    ### constructor
    def __init__(self,
                 sh    = shape(),
                 knots = np.zeros(shape=(0,1),dtype=object),
                 cp    = None ,
                 prt   =  False,
                 properties = {"periodic":None,"dtype":float}
                ):
        
        # making sure knots is a np.ndarray
        #knots = np.asarray([knots])
        # call init_knot
        self  = self.init_knot(sh,knots,periodic)
        
        # decido se stampare a schermo tutto 
        self._print = prt
        
        #tipo di variabile
        self.properties = properties
        #self.properties["dtype"] = dtype
           
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
    def init_knot(self,sh,knots,periodic):
        
        # assign value to class members
        #dim    = sh.dim()        
        self._sh = sh
        
        self._periodic = np.full(sh.dim(),False,dtype=bool)
        if periodic is None :
            periodic = self._periodic
        if self.dim() == 1:
            self._periodic[0] = periodic
        elif len(periodic) != sh.dim():
            print("Error : periodic array of wrong lenght")
            raise Exception()  
        else :
            self._periodic = periodic
        
        #check dimension
        if len(knots) > sh.dim() :
            print("Warning : knots lenght too long")
        elif len(knots) < sh.dim() :
            print("Error : knots lenght too short")
            raise Exception()          flo  
         
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
            self._cp = self._cp.astype(self.dtype)
        
        #
        self._stiffness_matrix = None
        self._overlap_matrix = None        
        #
        self._stiffness_matrix_BEM = None
        self._slp_matrix_BEM = None
        self._sol_BEM = None
        self._ind_sol_BEM = None
        self._lv_BEM = None        
        
        #
        self._ready_sm = False #stiffness matrix
        self._ready_om = False #overlap matrix        
        #
        self._ready_sm_BEM = False
        self._ready_slp_BEM = False
        self._ready_sol_BEM = False
        self._ready_lv_BEM = False
        self._ready_ind_sol_BEM = False
        
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
    
    ###
    def is_periodic(self):
        return np.any(self._periodic)
    
    ###
    def periodicity(self):
        index = self.control_points().index
        index = [ ii for ii in index]
        cp = pd.DataFrame(index = index,columns=["periodic"])
        cp["periodic"] = None
        for i in index:
            j = self.get_periodic_index(i)
            if j != i :
                cp.at[i,"periodic"] = j
        return cp
    ###
    def get_periodic_index(self,index):
        #https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-closed.html
        # 0 -> n-p
        # 1 -> n-p+1
        # ...
        # p-1 -> n-1
        # in general
        # first p
        # j -> t = n-p+j
        # last p
        # j -> j-n+p
        if self.dim() == 1 and type(index) is not tuple:
            index = [index]
        out = list(np.full(self.dim(),0))
        for i in range(self.dim()):
            j = index[i]
            n = self._kv[i].n()
            p = self._kv[i].p()
            if j <= p :
                out[i] = n-p+j-1#n-1-j
            elif j >= n-p-1:
                out[i] = j-n+p+1
            else :
                out[i] = j
        #if self.dim() != 1 :
        return tuple(out)
        #else :
        #    return out[0]
        
    ###
    def dof(self):
        il = self.index_list()
        it = [tuple(i) for i in il ]
        dof = self.periodicity()
        per = dof.copy()
        for i in it :
            j = per.at[i,"periodic"]
            #print(j)
            if j is not None and i in dof.index :
                dof.drop(j,inplace = True)
        return dof
    
    ### set control point value/coordinates
    def set_cp (self,index,value,check=True) :
        
        if self.codim() > 1 :
            value = np.asarray(value)
            value.reshape((self.codim(),1))
            
        if self.dim() != 1 :
            ti = tuple(index)
            self._cp[ti] = value
            if self.is_periodic():
                tip = self.get_periodic_index(index)
                self._cp[tip] = value
        else :
            self._cp[index] = value
            if self.is_periodic():
                tip = self.get_periodic_index(index)
                self._cp[tip] = value
                
        self._ready_sm = False
        self._ready_om = False
        self._ready_sm_BEM = False
        self._ready_slp_BEM = False
        self._ready_sol_BEM = False
        self._ready_lv_BEM = False
        self._ready_ind_sol_BEM = False        
        
        #self_ready_lv = False
        return True                                     
    ### some function returing some variables, not very usefull...
    def knots_vectors(self): 
        return self._kv.copy()
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
    ### da rivedere
    def clear_cp (self) : 
        cp = np.zeros(self._cp.shape,dtype=object)#,np.zeros(bs.codim()))
        cp.fill(np.zeros(self.codim()))
        if self.codim() == 1 :
            cp = cp.astype(float)
        self._cp = cp#np.full(self._cp.shape,np.full(self.codim(),0)).astype(float)
        self._ready_sm = False
        self._ready_om = False
        self._ready_sm_BEM = False
        self._ready_slp_BEM = False
        self._ready_sol_BEM = False
        self._ready_lv_BEM = False
        self._ready_ind_sol_BEM = False        
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
        curve      = Bspline(new_sh,curve_knot,dtype=self.dtype,periodic=self._periodic)
        
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
        return Bspline(out_sh,[sub_kv],surface_cp,dtype=self.dtype,periodic=self._periodic)    
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
    ### da rivedere
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
                    der_kv.append(knot_vector(kv.p(),kv.n(),kt,check_open=False)) 
                else :     
                    der_kv.append(knot_vector(kv.p()-1,kv.n()-1,kt[1:-1],check_open=False))         
            #der_kv = knot_vector(kv.p()-1,kv.n()-1,kt[1:-1]) # for i in self._kv]                
            derK = Bspline(der_sh,der_kv,periodic=self._periodic)

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
    
    ### da rivedere
    def basis_range(self):
            
        basis_range = list()#np.zeros(self.dim(),dtype=object)
        for k in range(0,self.dim()):
            kv = self._kv[k]
            kt = kv.knots()
            n  = kv.n()
            p  = kv.p() 
            data = {"index":range(0,n),                   "i-min":range(0,n),                   "i-max":np.arange(0,n)+p+1,                   "x-min":[ max(i,kt[p]) for i in kt[0:n]],                   "x-max":[ min(i,kt[n-p+1]) for i in kt[p+1:p+1+n]] }
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
    
    ### da rivedere
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
    
        #if self.is_periodic() == False:
        #    return am
        #else :
        #    index = [ list(i)[0] for i in am.index]
        #    #new_index = [ self.get_periodic_index(i) for i in index ]
        #    #ni = [ tuple([i]) for i in new_index]
        #    
        #    ni = [ self.get_periodic_index(i) for i in index ]
        #
        #    #adjacency matrix of periodic control points
        #    ams = am.copy()
        #    ams.index = ni
        #    ams.reindex(am.index,axis=0)
        #    out = ams+am
        #    out = out.astype(bool)
        #    return out
    
    ### da rivedere
    def basis_max_min(self,r,br=None):
        
        if br is None :
            br = self.basis_range()
        
        a = list()
        for k in range(0,self.dim()):
            a.append( [ br.at[r,("min",k)] , br.at[r,("max",k)] ] )    
            
        return a
    
    ### da rivedere
    def basis_overlap(self,r,c,br=None):
        
        if br is None :
            br = self.basis_range()
        
        rmm = self.basis_max_min(r,br)
        cmm = self.basis_max_min(c,br)
        o = list()
        for k in range(0,self.dim()):
            o.append( overlaps( rmm[k] ,cmm[k] ) )
            
        return o
    
    ### da rivedere
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

                # provo a modificare
                #X = [ np.delete(np.linspace(ov[k][0],ov[k][1],\
                #                            opts["delta"][k]+1,endpoint=False),0) \
                #     for k in range(0,self.dim()) ]
                X = np.full(self.dim(),0.0,dtype=object)
                for k in range(0,self.dim()):
                    # opts["delta"][k] = numero di punti interni
                    punti = np.linspace(ov[k][0],ov[k][1],opts["delta"][k]+1,endpoint=False)
                    X[k] = np.delete(punti,0)

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
        if "return_both" not in opts :
            opts["return_both"] = False
        #if "norm" not in opts :
        #    opts["norm"] = "L1"
            
        #if "del-edge" not in opts :
        #    opts["del-edge"] = True
        if opts2 is not None :
            opts.update(opts2)
        #if "ready_om" not in opts :
        #    opts["ready_om"] = True
        #if "ready_sm" not in opts :
        #    opts["ready_sm"] = True
            
        #if "ready_sm_BEM" not in opts :
        #    opts["ready_sm_BEM"] = True
        #if "ready_slp_BEM" not in opts :
        #    opts["ready_slp_BEM"] = True
        #if "ready_sol_BEM" not in opts :
        #    opts["ready_sol_BEM"] = True   
        #if "ready_ind_sol_BEM" not in opts :
        #    opts["ready_ind_sol_BEM"] = True   
        #if "ready_lv_BEM" not in opts :
        #    opts["ready_lv_BEM"] = True   
        
        varTrue = ["ready_sm_BEM","ready_slp_BEM","ready_sol_BEM",                   "ready_ind_sol_BEM","ready_lv_BEM",                   "copy_sm_BEM","copy_slp_BEM","copy_sol_BEM",                  "copy_ind_sol_BEM","copy_lv_BEM",                  "interpolation",                  "del-edge","ready_om","ready_sm"]
        
        for v in varTrue :            
            if v not in opts :
                opts[v] = True
            
            
        #if "interpolation" not in opts:
        #    opts["interpolation"] = True
        #if "ready_lv" not in opts :
        #    opts["ready_lv"] = False
        if "ready_trace" not in opts :
            opts["ready_trace"] = [ False for i in range(self.dim())]
        
        return opts
    
    ### da rivedere
    def edge(self):
        il = self.index_list()
        df = pd.DataFrame( il , index = tuple(il) ,columns = np.arange(0,self.dim()) )
        df["edge"] = True
        df["corner"] = True
        allx = np.arange(0,self.dim())
        e = np.zeros(self.dim())
        for i in range(0,df.shape[0]):
            for k in allx :                
                kv = self._kv[k]
                e[k] = df.iloc[i,k] <= 0 or df.iloc[i,k] >= kv.n()-1
                #e[k] = df.iloc[i,k] < kv.p() or df.iloc[i,k] > kv.n()-kv.p()-1
            df.at[df.index[i],"edge"] = np.any(e)
            df.at[df.index[i],"corner"] = np.all(e)
        df.drop(columns=allx,inplace=True)
        return df       
    
    ### da rivedere
    def approximate(self,func,opts=None):
        #http://hplgit.github.io/INF5620/doc/pub/sphinx-fem/._main_fem002.html
        
        opts = self.prepare_opts(opts)

        #
        index = self.index_list()
        it = [ tuple(j) for j in index ]        

        #edge -> da rivedere
        edge = self.edge()

        #overlap matrix
        om = self.overlap_matrix(opts)
        om.replace(np.nan,0.0,inplace=True)

        if self.dim() > 1 :

            #indici dei dof interni e di bordo
            index_int  = edge.index[ edge["edge"] == False ]
            index_edge = edge.index[ edge["edge"] == True  ]

            #degrees of freedom: internal
            dof_int = om.copy()
            dof_int.drop( index_edge ,inplace = True,axis=0)
            dof_int.drop( index_edge ,inplace = True,axis=1)

            #degrees of freedom: edge
            dof_edge = om.copy()
            dof_edge.drop( index_edge ,inplace = True,axis=0)
            dof_edge.drop( index_int  ,inplace = True,axis=1)

            #load vector
            lv = self.load_vector(func,opts)        
            lv.replace(np.nan,0.0,inplace=True)
            lv.drop( index_edge  ,inplace = True)  

            #convert into numpy array
            dinp = np.asarray(dof_int)
            denp = np.asarray(dof_edge)            
            
            #
            gDv = self.Dirichlet_BC(func,opts)
                        
            # -> da rivedere
            # qui devo introdurre del codice per gestire 
            # il fatto che la Bspline può essere periodica

            #load vector
            lvnp = np.asarray(lv)
            #punti di bordo
            gDnp = np.asarray(gDv,dtype=self.dtype)
            #prodotto righe per colonne
            edlv = np.dot(denp,gDnp)

            #solve linear system
            cpint = np.linalg.solve(dinp,lvnp-edlv) 

            #preparo gli indice della variabile di output
            #index = self.index_list()
            #it = [ tuple(j) for j in index ]
            #out = pd.DataFrame(index=it,columns=["cp"])

            # -> da rivedere
            # cpint.index sarà diverso da index_int
            # sostituire index_int con cpint.index
            
            #assegno ai control points i valori calcolati
            #valori interni
            for i in range(len(index_int)):
                j = index_int[i]
                self._cp[j] = cpint[i]
                #out.iloc[out.index == j] = cpint[i]
                
            # -> da rivedere
            # gDnp.index sarà diverso da index_edge
            # sostituire index_edge con gDnp.index

            #assegno ai control points i valori calcolati
            #valori di bordo interpolanti
            for i in range(len(index_edge)):
                j = index_edge[i]
                self._cp[j] = gDnp[i]
                #out.iloc[out.index == j] = gDnp[i]


        else :                

            if opts["interpolation"] == True :
            #indici dei dof interni e di bordo
                index_int  = edge.index[ edge["corner"] == False ]
                index_edge = edge.index[ edge["corner"] == True  ]
            else :
                index_int  = edge.index#[ edge["corner"] == False ]
                index_edge = []#edge.index[ edge["corner"] == True  ]

            #degrees of freedem: internal
            dof_int = om.copy()
            dof_int.drop( index_edge ,inplace = True,axis=0)
            dof_int.drop( index_edge ,inplace = True,axis=1)

            #degrees of freedem: edge            
            dof_edge = om.copy()
            dof_edge.drop( index_edge ,inplace = True,axis=0)
            dof_edge.drop( index_int  ,inplace = True,axis=1)

            #convert into numpy array
            dinp = np.asarray(dof_int)
            denp = np.asarray(dof_edge)

            #load vector: dof interni
            #print("load vector")
            lv = self.load_vector(func,opts)
            lv.drop( index_edge  ,inplace = True) 


            #
            XY = pd.DataFrame(index = index_edge,columns=np.arange(0,self.dim()))
            for k in range(self.dim()):
                kv = self.knots_vectors()[k]
                nmax = kv.n()-1
                xmin = min(kv.knots())
                xmax = max(kv.knots())
                for i in index_edge:
                    if i[k] == nmax :
                        XY.at[i,k] = xmax
                    else :
                        XY.at[i,k] = xmin

            xy = np.asarray(XY).astype(float)
            #if self.dim() == 1 :
            #    xy = xy.reshape((len(xy),))

            #dof di bordo
            #print("gDnp")
            gDnp = func(xy)#.astype(float)

            ###        
            #print("ciao")
            lvnpND = np.asarray(lv["cp"])#.astype(float)
            lvnpND = np.zeros(shape=(len(lv),self.codim()))
            for i in range(len(lv)):
                lvnpND[i,:] = lv["cp"][i]    

            if self.codim() == 1 :
                a = gDnp.copy()
                gDnp = np.zeros(shape=(len(a),self.codim()))
                gDnp[:,0] = a
                del a
            gDnpND = gDnp

            #prodotto righe per colonne
            #edlv = np.dot(denp,gDnp)
            #edlvND = edlv
            
            # -> da rivedere
            # qui devo introdurre del codice per gestire 
            # il fatto che la Bspline può essere periodica


            out = pd.DataFrame(index=index_int,columns=np.arange(0,self.codim()))
            #index = self.index_list()
            #it = [ tuple(j) for j in index ]
            for k in range(self.codim()):

                gDnp = gDnpND[:,k]
                edlv = np.dot(denp,gDnp)
                #edlv = edlvND[:,k]
                lvnp = lvnpND[:,k]


                gDnp = gDnp.reshape((len(gDnp),))
                edlv = edlv.reshape((len(edlv),))

                #
                cpint = np.linalg.solve(dinp,lvnp-edlv) 

                #cp = pd.DataFrame(index=it,columns=["cp"])

                out[k] = cpint

            #assegno ai control points i valori calcolati
            #valori interni
            for i in range(len(index_int)):
                j = index_int[i]
                self._cp[j] = np.asarray(out.iloc[i,:])
                #out.iloc[out.index == j] = cpint[i]

            #assegno ai control points i valori calcolati
            #valori di bordo interpolanti
            for i in range(len(index_edge)):
                j = index_edge[i]
                self._cp[j] = gDnpND[i]
                #out.iloc[out.index == j] = gDnp[i]
        
        return self.control_points()
    
    ##
    def _scalar(self):#restituisce la stessa Bspline ma con cp scalari
        sh = shape(dim=self.dim(),codim=1)
        kv = self._kv        
        return Bspline(sh,kv,periodic=self._periodic)#cp: default value
    
    ###
    def get_gD_shape(self):
        index = self.index_list("edge")
        it = [ tuple(j) for j in index ]
        columns = ["Dirichlet"]
        return pd.DataFrame(data=np.full(len(it),0),index=it,columns=columns)
        
    ### da rivedere: controllare
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
        kv     = [ self.knots_vectors()[ii] for ii in range(0,self.dim()) if ii != n ]  
        per    = self._periodic.copy().drop(n)
        curve  = Bspline(new_sh,kv,periodic=per)
        
        self._ready_trace[n] = True
        self._trace_Bspline[n] = curve.copy()
        
        return curve
            
    ### da rivedere
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
       
    ### da rivedere
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
        gDnp = np.asarray(gDv,dtype=self.dtype)
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
            self._cp = self._cp.astype(self.dtype)
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
            #if opts["norm"] == "L1" :
            #modifica
            res = np.mean(y, axis=0) if self.codim() > 1 else np.mean(y)
            #elif opts["norm"] == "L2" :
            #modifica
            #    res = np.sqrt(np.sum(np.power(y,2.0),axis=0))/len(y) if self.codim() > 1 \
            #    else np.sqrt(np.sum(np.power(y,2.0)))/len(y)
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
    #BEM
    
    ##
    def stiffness_matrix_BEM(self,k=None,opts=None):
        
        opts = self.prepare_opts(opts)
        
        if opts["ready_sm_BEM"] == True and self._ready_sm_BEM == True :
            return self._stiffness_matrix_BEM

        #if opts["ready_sm_BEM"] == False:
        #    self._ready_sm_BEM = False

        #if self._ready_sm_BEM == True :
        #    return self._stiffness_matrix_BEM
        
        if opts["print"] :
            print("stiffness_matrix_BEM")

        conta = 1
        for k in range(0,self.dim()):
            conta = conta * opts["delta"][k]
        Xintegration = np.zeros(shape=(conta,self.dim()))
        Yintegration = Xintegration.copy()
        del conta

        #X = np.full(self.dim(),0.0,dtype=object)
        X = np.full(self.dim(),0.0,dtype=object)
        for k in range(0,self.dim()):
            X[k] =  np.full(opts["delta"][k],0.0)
        Y = X.copy()        

        #a djacency matrix
        #am = self.adjacency_matrix()
        il = self.index_list()
        it = [ tuple(i) for i in il ]
        # basis range
        br = self.basis_range()
        # derivative
        der = self.derivative()
        # basis function: left 
        scalar = self._scalar()
        # basis function: right 
        left   = scalar.copy()
        left.clear_cp()
        #
        right  = scalar.copy()        
        right.clear_cp()
        #
        del scalar

        # out variable
        out = pd.DataFrame(0.0,index=it,columns=it,dtype=np.complex)

        # unità immaginaria
        I = np.complex(0,1)
        
        def foundamental(d):
            return scipy.special.hankel1(np.full(len(d),0),k*d)*I/4.0  

        def integrate(xy):
            x = xy[:,0]
            y = xy[:,1]
            l = left.evaluate(x)
            r = right.evaluate(y)#.conjugate()
            dl = [ norm(i) for i in der.evaluate(x) ]
            dr = [ norm(i) for i in der.evaluate(y) ]
            xx = self.evaluate(x).astype(self.dtype)
            yy = self.evaluate(y).astype(self.dtype)
            d = np.asarray([ norm(i) for i in xx-yy ])
            f = foundamental(d)
            return f*d*l*dl*dr


        ###
        # CICLO
        successful = False
        res =  np.complex(0,0)
        n = len(it)#am.shape[0]
        conta = 0
        tot = int(n*(n+1)/2)
        for i in range(0,n):

            #
            r = it[i]#am.index[i] 
            left.set_cp(r,1.0)
            der_left  = left.derivative()

            areaX = 1.0
            for k in range(0,self.dim()):
                # ATTENZIONE ALLA DERIVATA LOGARITMICA
                xmin = br.at[r,("min",k)]
                xmax = br.at[r,("max",k)]
                delta = xmax-xmin
                punti0 = np.linspace(xmin,xmax,opts["delta"][k]+1,endpoint=True)
                puntiX  = np.delete(punti0,0)
                #numradX = (0.5-random.rand(len(punti)))*delta/(opts["delta"][k]+2)
                #numradY = (0.5-random.rand(len(punti)))*(ov[k][1]-ov[k][0])/(opts["delta"][k]+2)
                #X[k] = punti+numradX
                areaX = areaX*delta
                #Y[k] = punti+numradY

            for j in range(i,n):

                c = it[j]#am.columns[j]
                
                if opts["print"] :
                    print(conta,"/",tot,end="\r")

                # LA MATRICE E' DENSA
                #if am.at[r,c] is False :
                #    continue

                #
                right.set_cp(c,1.0)
                #ov = self.basis_overlap(r,c,br)
                der_right = right.derivative()

                # ATTENZIONE:
                # preparare X e Y in due step diversi, si riferiscono a due funzioni di base diverse
                # una dentro e l'altra fuori dal ciclo
                # considerare i range di entrambe le funzioni di base
                # provare a calcolare X e Y all'inizio di tutto quanto
                
                
                # questi sono da modificare  
                areaY = 1.0
                for k in range(0,self.dim()):
                    # ATTENZIONE ALLA DERIVATA LOGARITMICA
                    xmin = br.at[c,("min",k)]
                    xmax = br.at[c,("max",k)]
                    delta = xmax-xmin
                    punti0 = np.linspace(xmin,xmax,opts["delta"][k]+1,endpoint=True)
                    puntiY  = np.delete(punti0,0)
                    areaY = areaY*delta
                 
                
                successful = False
                while not successful :
                    
                    for k in range(0,self.dim()):
                        # ATTENZIONE ALLA DERIVATA LOGARITMICA
                        numradX = (0.5-random.rand(len(puntiX)))*delta/(opts["delta"][k]+2)
                        numradY = (0.5-random.rand(len(puntiY)))*delta/(opts["delta"][k]+2)
                        X[k] = puntiX+numradX
                        Y[k] = puntiY+numradY                     
                    

                    # devo usare np.meshgrid
                    mX = np.meshgrid(*X)
                    mY = np.meshgrid(*Y)
                    for k in range(0,self.dim()):
                        Xintegration[:,k] = mX[k].flatten()
                        Yintegration[:,k] = mY[k].flatten()

                    # questo è ottimizzabile
                    meshN = len(Xintegration)*len(Yintegration)
                    mesh = np.zeros((meshN,2),dtype=object)
                    mesh.fill(Xintegration[0])
                    n1 = len(Xintegration)
                    n2 = len(Yintegration)
                    for i in range(n1):
                        for j in range(n2):
                            k = n1*i+j
                            mesh[k,0] = Xintegration[i]
                            mesh[k,1] = Yintegration[j]
                    # piccola modifica
                    #mesh = mesh [ mesh[:,0] != mesh[:,1]  ]

                    # piccola modifica
                    if self.dim() == 1:
                        mesh = mesh.astype(float)

                    ###
                    # VALUTAZIONE DELL'INTEGRALE

                    y = integrate(mesh)
                    res = np.mean(y)*areaX*areaY
                    
                    #
                    successful = not np.isnan(res)
                    
                    if not successful :
                        print("found nan value for [",r,",",c,"]: repeating the cycle")

                    ###
                    # FINE
                    
                   
                #    
                out.at[r,c] = res
                out.at[c,r] = res #matrice simmetrica

                #cancello
                right.set_cp(c,0.0)
                
                #incremento
                conta = conta +1

                #cancello
            left.set_cp(r,0.0)
         
        # TENGO CONTO DELLA PERIODICITA'
        # le funzioni di base non corrispondono ai control points
        #
        #il = self.index_list()
        #it = [tuple(i) for i in il ]
        dof = self.dof()
        
        #self.periodicity()
        #per = dof.copy()
        #for i in it :
        #    j = per.at[i,"periodic"]
        #    #print(j)
        #    if j is not None and i in dof.index :
        #        dof.drop(j,inplace = True)


        persm = out.copy()
        delindex = [ i for i in out.index if i not in dof.index ]
        perindex = [ self.get_periodic_index(i) for i in delindex  ]
        index = out.index
        for i in index:
            for j in perindex :
                #ii = self.get_periodic_index(i)
                jp = self.get_periodic_index(j)
                persm.at[i,j] += persm.at[i,jp]

        for i in index:
            for j in perindex :
                #ii = self.get_periodic_index(i)
                jp = self.get_periodic_index(j)
                persm.at[j,i] += persm.at[jp,i]
                #persm.at[j,i] = persm.at[i,j]

        for i in delindex:
            persm.drop(i,inplace=True,axis=0)
            persm.drop(i,inplace=True,axis=1)
                    
        #    
        if opts["copy_sm_BEM"] == True :
            self._stiffness_matrix_BEM = persm.copy()
            self._ready_sm_BEM = True
        if opts["print"] :
            print("Finished")
        if opts["return_both"]:
            return persm,out
        return persm,out
    
    ###
    def load_vector_BEM(self,gD=None,opts=None):
        
        opts = self.prepare_opts(opts)
        
        if opts["ready_lv_BEM"] == True and self._ready_lv_BEM == True :
            return self._lv_BEM
            
        #if opts["ready_lv_BEM"] == False:
        #    self._ready_lv_BEM = False

        #if self._ready_lv_BEM == True :
        #    return self._lv_BEM
        
        #edge = self.edge()

        #
        X = np.full(self.dim(),0.0,dtype=object)
        for k in range(0,self.dim()):
            X[k] =  np.full(opts["delta"][k],0.0)
        conta = 1
        for k in range(0,self.dim()):
            conta = conta * opts["delta"][k]
        Xintegration = np.zeros(shape=(conta,self.dim()))

        #
        br = self.basis_range()
        il = self.index_list()
        il = [ tuple(i) for i in il ]
        der = self.derivative()
        scalar = self._scalar()
        scalar.clear_cp()
        out = pd.DataFrame(index = il, columns=["cp"],dtype=object)
        #out["index"] = il
        #out.set_index("index",inplace=True)

        #definisco la norma di Lebesgue

        def integrate(x):
            a = scalar.evaluate(x)
            xx = self.evaluate(x)
            b = gD(xx)
            d0 = der.evaluate(x)
            d = [ norm(i) for i in d0 ]
            return a*b*d

        for i in il :

            #i = tuple(i)
            scalar.set_cp(i,1.0)
            #ov = self.basis_overlap(i,i,br) #overlap

            areaX = 1.0
            for k in range(0,self.dim()):
                xmin = br.at[i,("min",k)]
                xmax = br.at[i,("max",k)]
                delta = xmax-xmin
                punti0 = np.linspace(xmin,xmax,opts["delta"][k]+1,endpoint=True)
                punti  = np.delete(punti0,0)
                numradX = (0.5-random.rand(len(punti)))*delta/(opts["delta"][k]+2)
                X[k] = punti+numradX
                areaX = areaX*delta

            m = np.meshgrid(*X)
            for k in range(0,self.dim()):
                Xintegration[:,k] = m[k].flatten()

            y = integrate (Xintegration)
            out.at[i,"cp"] = np.mean(y) * areaX

            scalar.set_cp(i,0.0)
            
        # TENGO CONTO DELLA PERIODICITA'
        dof = self.dof()
        persm = out.copy()
        delindex = [ i for i in out.index if i not in dof.index ]
        perindex = [ self.get_periodic_index(i) for i in delindex  ]
        index = out.index
        for i in index:
            ip = self.get_periodic_index(i)
            persm.at[i,"cp"] += persm.at[ip,"cp"]
        for i in delindex:
            persm.drop(i,inplace=True,axis=0)

            
        if opts["copy_lv_BEM"] == True :
            self._lv_BEM = persm.copy()
            self._ready_lv_BEM = True
        if opts["return_both"]:
            return persm,out
        return persm
    
    ###
    def single_layer_potential_basis_BEM(self,XY=None,k=None,opts=None):
        # variabile di output
        # matrice con:
        # - righe: punti x dove valutare la soluzione
        # - colonne : funzioni di base
        #out = pd.DataFrame(index=np.arange(0,len(XY)),columns=lv.index)
        opts = self.prepare_opts(opts)

        if opts["ready_slp_BEM"] == False:
            self._ready_slp_BEM = False

        if self._ready_slp_BEM == True :
            return self._slp_matrix_BEM
        
        #il = self.index_list()
        #il = [ tuple(i) for i in il ]
        dof = self.dof()
        outnp = np.zeros(shape=(len(XY),len(dof)),dtype=object)

        #giusto per definirlo
        x0 = XY[0,:]
        # foundamental solution with x fixed
        I = np.complex(0,1)
        def foundamental_x(y):
            d = np.asarray([ norm(i)  for i in y-x0 ])
            return scipy.special.hankel1(np.full(len(d),0),k*d)*I/4.0 

        opts2 = opts
        opts2["ready_lv_BEM"] = False
        opts2["copy_lv_BEM"] = False
        lenXY = len(XY)
        for i in range(lenXY):
            if opts["print"] :
                print(i,"/",lenXY,end="\r")
            x0 = XY[i,:]
            lvx0 = self.load_vector_BEM(foundamental_x,opts2)
            outnp[i] = lvx0["cp"]
        if opts["print"] :
            print("Finished")

        out = pd.DataFrame(data=outnp,index = [ tuple(i) for i in XY] ,columns=dof.index)#["x",il])
        
        if opts["copy_slp_BEM"] == True :
            self._slp_matrix_BEM = out.copy()
            self._ready_slp_BEM = True
        return out
    
    ###
    def single_layer_potential_BEM(self,sol,slpB,XY,k,opts=None):
        
        # sol  : solution of linear system
        # slpB : single layer potential (basis)
        
        opts = self.prepare_opts(opts)

        if opts["ready_sol_BEM"] == False:
            self._ready_sol_BEM = False

        if self._ready_sol_BEM == True :
            return self._sol_BEM
        
        #        
        slpBnp = np.asarray(slpB)        
        solnp = np.asarray(sol["value"])

        #
        outnp = np.dot(slpBnp,solnp)

        #
        out = pd.DataFrame(index=slpB.index,columns=["x","value"])
        out["x"] = [ tuple(i) for i in XY]
        out["value"] = outnp
        
        if opts["copy_sol_BEM"] == True :
            self._sol_BEM = out.copy()
            self._ready_sol_BEM = True
        return out
    
    ###
    def indirect_solution_BEM(self,sm=None,lv=None,opts=None):
        
        opts = self.prepare_opts(opts)

        if opts["ready_ind_sol_BEM"] == False:
            self._ready_ind_sol_BEM = False

        if self._ready_ind_sol_BEM == True :
            return self._ind_sol_BEM
        
        # convert into numpy array
        smnp = np.asarray(sm).astype(np.complex)
        lvnp = np.asarray(lv).reshape((len(lv),)).astype(np.complex)
        
        #
        sol = scipy.linalg.solve(smnp,lvnp,assume_a="sym") #lo salvo su file
        
        out = pd.DataFrame(data=sol,index=sm.index,columns=["value"])
     
        self._ind_sol_BEM = out.copy()
        self._ready_ind_sol_BEM = True
        return out
    
    ###
    def internal_points(self,XY,NN,xmin,xmax,opts):
        
        t = np.linspace(xmin,xmax,NN,endpoint=True)
        xy   = self.evaluate(t)
        area = max(t)-min(t)

        def to_complex(xy):
            return np.asarray([ np.complex(i[0],i[1]) for i in xy ])

        z = to_complex(xy)
        der = self.derivative()
        d = der.evaluate(t)
        jac = to_complex(d)#np.asarray([np.sqrt(np.sum(np.power(i,2.0))) for i in d ])

        I = np.complex(0,1)

        x0 = to_complex(XY)
        wn = np.zeros(x0.shape,dtype=np.complex)

        for i in range(len(x0)):
            delta = z-x0[i]
            integrand = jac /delta * area
            wn[i] = np.mean(integrand)/(2*np.pi*I)
            
        absolute = np.absolute(wn)
        out = absolute > 0.5
        return out
        
    ###
    def BEM(self,uinc=None,k=None,XY=None,opts=None):
        
        opts = self.prepare_opts(opts)

        #indirect solution: Galerkin  method
        sm = self.stiffness_matrix_BEM(k=k,opts=opts) #lo salvo su file
        lv = self.load_vector_BEM(gD=uinc,opts=opts) #lo salvo su file        
        sol = self.indirect_solution_BEM(sm,lv,opts)
        
        #
        slpB = self.single_layer_potential_basis_BEM(XY=XY,k=k,opts=opts) #lo salvo su file        
                      
        #
        slp = self.single_layer_potential_BEM(sol=sol,slpB=slpB,XY=XY,k=k,opts=opts)
                
        # convert into numpy array
        Xnp,Valnp = self.sol_to_np_BEM(slp)
        
        #self._sol_BEM = slp.copy()
        #self._ready_sol_BEM = True
        return slp,Xnp,Valnp
    
    ###
    def sol_to_np_BEM(self,sol):
        # convert into numpy array
        r = len(sol)
        c = len(sol["x"][0])

        Xnp = np.zeros(shape=(r,c))

        for i in range(r):
            for j in range(c):
                Xnp[i,j] = sol["x"][i][j]

        #
        Valnp = np.asarray(sol["value"],dtype=np.complex)
        return Xnp,Valnp
                    
    ###
    def save(self,variable,filename):
        if variable == "sm":
            var = self._stiffness_matrix
        elif variable == "om":
            var = self._overlap_matrix
        elif variable == "lv":
            var = self._load_vector
        elif variable == "cp" :
            var = self.control_points()
        elif variable == "sm-BEM" :
            var = self._stiffness_matrix_BEM
        elif variable == "slp-BEM" :
            var = self._slp_matrix_BEM
        elif variable == "lv-BEM" :
            var = self._lv_BEM
        elif variable == "ind_sol-BEM" :
            var = self._ind_sol_BEM  
        elif variable == "sol-BEM" :
            var = self._sol_BEM  
        var.to_csv(filename,index_label="index")
        
    ###
    def load(self,variable,filename):
        
        var = pd.read_csv(filename)
                
        #l eggo e preparo gli indici
        int_index = ["sm","sm-BEM","om","cp","ind_sol-BEM","lv-BEM"]
        if variable in int_index:
            var.index = [tuple_from_str(i) for i in var["index"] ]
            var = var.drop('index',axis=1)
            
                
        # leggo e preparo le colonne
        int_col = ["sm","sm-BEM","om"]
        if variable in int_col:
            var.columns = [tuple_from_str(i) for i in var.columns]
        
        
        # leggo e preparo i valori
        if variable == "sm":  
            self._stiffness_matrix = var.copy()
            self._ready_sm = True     
                
        elif variable == "sm-BEM":
            var = var.applymap(np.complex)
            self._stiffness_matrix_BEM = var.copy()
            self._ready_sm_BEM = True 
                           
        elif variable == "om":
            self._overlap_matrix = var.copy()
            self._ready_om = True     
            
        elif variable == "slp-BEM":
            var.index = [ tuple_float_from_str(i) for i in var["index"]]
            var = var.drop('index',axis=1)
            var.columns = [tuple_from_str(i) for i in var.columns]
            var = var.applymap(np.complex)
            self._slp_matrix_BEM = var.copy()
            self._ready_slp_BEM = True   
        
        elif variable == "sol-BEM":
            var.index = [ tuple_float_from_str(i) for i in var["index"]]
            var["x"] = [ tuple_float_from_str(i) for i in var["x"]]
            var = var.drop('index',axis=1)
            var["value"] = [np.complex(i) for i in np.asarray(var["value"])]
            self._sol_BEM = var.copy()
            self._ready_sol_BEM = True   
            Xnp,Valnp = self.sol_to_np_BEM(var)
            return var,Xnp,Valnp
        
        elif variable == "ind_sol-BEM":
            var["value"] = [np.complex(i) for i in np.asarray(var["value"])]
            self._ind_sol_BEM = var.copy()
            self._ready_ind_sol_BEM = True     
            
        elif variable == "lv-BEM":
            var["cp"] = [np.complex(i) for i in np.asarray(var["cp"])]
            self._lv_BEM = var.copy()
            self._ready_lv_BEM = True     
                     
        elif variable == "cp" :
            var.columns = [ int(i) for i in var.columns ]
            index = var.index
            for i in index:
                for k in range(self.codim()):
                    self._cp[i,k] = var.at[i,k]
            #self._cp = var.copy()
        
        return var

###
def norm(x):
    return np.sqrt(np.sum(np.power(x,2.0)))        
      
###
def tuple_from_str(string,comand=r'[0-9]+'):
    return tuple(map(int, re.findall(comand, string)))

###
def tuple_float_from_str(string,comand=r"[+-]?\d+(?:\.\d+)?\.[+-]?\d+(?:\.\d+)?"):
    return tuple(map(float, re.findall(comand, string)))

    
###        
def overlaps(a, b):
    c = [ max(a[0],b[0]) , min(a[1],b[1]) ]
    
    if c[0] >= c[1] :
        return None
    else :
        return c


# 
# 
