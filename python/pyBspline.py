#!/usr/bin/env python
# coding: utf-8

# In[10]:


get_ipython().system('jupyter nbconvert --to script pyBspline.ipynb')


# ## Shape class

# In[31]:


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

# In[32]:


import numpy as np

class knot_vector ():
    
    #_pol_order  = 0 #grado polinomiale
    #_basis_card = 0 #cardinalità della base 
    #_vect       = np.ndarray(shape=(0,))
      
    ###
    def __init__(self,
                 p=0,
                 n=0,
                 cls=np.zeros(shape=(0,1))):
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
    
    ###
    def p     (self): return self._pol_order
    def n     (self): return self._basis_card
    def knots (self): return self._vect 
    ###
    def __len__(self): return len(self.knots())
    ###
    def show  (self): 
        print("polinomial degree : " , self.p())
        print("base caridnality  : " , self.n())
        print("knots             : " , self.knots())
    
                            


# ## Bspline class

# In[33]:


import copy
import pandas as pd
from scipy import integrate
import itertools 
import time

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
        return True                                     
    ### some function returing some variables, not very usefull...
    def knots_vectors(self): return self._kv
    ###
    def get_knots_vector (self,index=None) :
        if index == None :        
            return self._kv
        else :
            return self._kv[index]
    ###
    def control_points (self) : 
        return self._cp    
    ###
    def clear_cp (self) : 
        self._cp = np.zeros(self._cp.shape)
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
            if t[i] <= x and x < t[i+1] :
                output = i
        #if x == t[-1] : #ultimo elemento
        #    output = i
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
        c = self.control_points()
        
        #
        k = self._find_k(x,t)                
        #print("k:",k)
        if k < 0 :
            if self._print == True :
                print("deBoor error : k<0")
            return self.Type_out()
        
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
        
        
        #print("k:",k)
        #print("x:",x)
        #print("t:",t)
        #print("c:",c)
        #print("type(c)",type(c))
        #print("p:",p)
        #print("len(c):",len(c))
        #print("len(t):",len(t))  
        
        
        #d = [ c[j + k - p] if for j in range(0, p+1)]
        d = np.zeros(p+1,dtype=object)
        #d.fill(self.Type_out())
        for  j in range(0,p+1) :
            d[j] = c[j+k-p] if 0 <= j+k-p < len(c) else self.Type_out()
            #if 0 <= j+k-p < len(c) :
            #    d[j]=c[j+k-p]
            #else :
            #    d[j] = self.Type_out()        
        
        #print("d:",d)
        #print("type(d)",type(d))

        for r in range(1, p + 1):
            #print("r:",r)
            for j in range(p, r - 1, -1): 
                right = j+1+k-r
                left  = j+k-p
                if 0 <= left < len(t) and 0 <= right < len(t) :
                    #print("j:",j)
                    #print("j+k-p:",j + k - p)
                    #print("j+1+k-r:",j + 1 + k - r)
                    alpha = (x - t[left]) / (t[right] - t[left])
                    #print("alpha(",j,"):",alpha)
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
    def derivative(self,n=1):        
        
        if self.dim() == 1 :        
            #http://public.vrac.iastate.edu/~oliver/courses/me625/week5b.pdf        
            # n = ordine della derivata
            # calcolo la derivata n volte
            for j in range(0,n):     
                
                der_kv = list()
                for k in range(0,self.dim()):
                
                    der_sh = self._sh                
                    kv     = self._kv[k]
                    kt     = kv.knots()
                    p      = kv.p()
                    der_kv.append(knot_vector(kv.p()-1,kv.n()-1,kt[1:-1]))         
                #der_kv = knot_vector(kv.p()-1,kv.n()-1,kt[1:-1]) # for i in self._kv]                
                der    = Bspline(der_sh,der_kv)

                out = list()
                N = list()
                X = list()
                for k in range(0,self.dim()):
                    #prendo solo quelli lungo un asse
                    kv     = self._kv[k]
                    kt     = kv.knots()
                    p      = kv.p()
                    #devo ciclare su tutte le altre dimensioni
                    N.append(kv.n())
                    X.append(kt)
                    
                    for i in range(0,kv.n()-1):
                        #print("i:",i)
                        #print("self.get_cp(i+1):",self.get_cp(i+1))
                        #print("self.get_cp(i)  :",self.get_cp(i))
                        #print("kt[i+2]         :",kt [i+2])
                        #print("kt[i+1]         :",kt[i+1])
                        cp = p*(self.get_cp(i+1) - self.get_cp(i)) / ( kt[i+p+1] - kt[i+1] ) 
                        der.set_cp(i,[cp]) 

                #if j != n-1 :
                self = der.copy()

            return der    
        else :
            
            
            #der.clear_cp()

            out = list()
            #derK = der.copy()
            for K in range(0,self.dim()):
                
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
                        derK.set_cp(ii,cp) 
                out.append(derK)

            return out

    ###
    def jacobian(self,x):
        #sistemo dimensioni di x
        X = self._correct_type_and_shape(x)
        #
        if len(X) == 1 :
            
            if self.dim() == 1 : #calcolo del gradiente
                
                der = self.derivative()
                return der.evaluate(X)
            
                #return self._deBoor(X,der=True)
            else : #calcolo della matrice jacobiana

                output = np.zeros(shape=(self.dim(),self.codim()))
                for i in range(self.dim()) :
                    newBspline = self._transpose(0,i)
                    curve_final = newBspline._iterative_deBoor(x)
                    output[i] = curve_final.jacobian(X)
        else :
            return self._correct_type_and_shape([ self.jacobian(j) for j in X ])
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
    def index_list(self):
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
        return index
    
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

        #adjacency matrix delle derivate
        am = self.adjacency_matrix()

        smd1D = am.copy()    
        smd1D[ smd1D == False ] = 0.0 #np.nan
        smd1D[ smd1D == True ] = 0.0 #np.nan

        #calcolo il prodotto scalare dei gradiente
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
                if opts["norm"] == "L1" :
                    res = np.mean(y)*area
                elif opts["norm"] == "L2" :
                    print("L2 norm")
                    res = None #np.sqrt(np.mean(y))*area

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
        return opts
    
    ###
    def edge(self):
        il = self.index_list()
        df = pd.DataFrame( il , index = tuple(il) ,columns = np.arange(0,self.dim()) )
        df["edge"] = True
        allx = np.arange(0,self.dim())
        e = np.zeros(self.dim())
        for i in range(0,df.shape[0]):
            for k in allx :                
                kv = self._kv[k]
                e[k] = df.iloc[i,k] < kv.p() or df.iloc[i,k] > kv.n()-kv.p()-1
            df.at[df.index[i],"edge"] = np.any(e)
        df.drop(columns=allx,inplace=True)
        return df
    
    ###
    def stiffness_matrix(self,opts=None):
        
        opts = self.prepare_opts(opts)

        #matrice di overlap delle derivate
        omd = list()
        der = self.derivative() if self.dim() > 1 else [self.derivative()]
        if opts["print"] == True : print("preparation",end="\r")
        for k in range(0,self.dim()):
            if opts["print"] == True : print("\ndimension :",k)  
            omd1D = der[k].overlap_matrix(opts)
            omd.append(omd1D)
        del omd1D
        #smd : overlap matrix derivatives
        #smd : matrice di overlap delle derivate

        ###
        #calcolo la stifness matrix vera e propria
        am = self.adjacency_matrix()
        #matrice di stiffness di una sola dimensione
        sm1D = pd.DataFrame(0.0,index=am.index,columns=am.columns)
        n = sm1D.shape[0]#quadrata
        
        # escludo subito dal calcolo i termini di bordo 
        #impostando a nan il valore in sm1D
        edge = self.edge()
                
        
        #sm1D = am.copy()    
        #sm1D[ sm1D == False ] = 0.0 #None

        #stiffness matrix: sum over all dimensione
        #matrice stiffness finale, è la somma della matrici parziali (su una sola dimensione)
        #questa è in realtà una lista contenente le matrici parziali
        sm = list()

        

        #creo una copia della Bspline
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
                if edge.at[r,"edge"] == True :
                    sm1D.at[r,:] = np.nan
                    sm1D.at[:,r] = np.nan
                    continue
                    
                #left.clear_cp()
                left.set_cp(r,1.0)
                #calcolo la derivata della funzione di base di sinistra
                dl = left.derivative()[k] if self.dim() > 1 else left.derivative()

                #cancello, tanto non mi serve valutarla
                #mi servono soltanto i coefficienti della derivata
                left.set_cp(r,0.0)


                #control points, modifico il tipo e cerco quelli non nulli
                cpl = dl._cp.astype(float) 
                #non zero control points (left) indices
                nzcpl = np.argwhere( cpl != 0.).tolist() 


                #ciclo su tutte le altre funzioni di base
                for j in range(i,n):
                    #if opts["print"] == True : print(i,"-",j,end="\r")

                    #indice della funzione di base di destra (c=colonne)
                    c = am.columns[j]   
                    
                    #controllo che il termine non sia di bordo
                    if edge.at[c,"edge"] == True :
                        sm1D.at[c,:] = np.nan
                        sm1D.at[:,c] = np.nan
                        continue
                    
                    
                    #right.clear_cp()
                    right.set_cp(c,1.0)
                    dr = right.derivative()[k] if self.dim() > 1 else right.derivative()
                    #cancello, tanto non mi serve valutarla
                    #mi servono soltanto i coefficienti della derivata
                    right.set_cp(c,0.0)

                    #control points, modifico il tipo e cerco quelli non nulli
                    cpr = dr._cp.astype(float)
                    #non zero control points (right) indices
                    nzcpr = np.argwhere( cpr != 0.).tolist() 

                    #attenzione all'ordine
                    #genero tutte le coppie
                    allpairs = list(itertools.product(nzcpl,nzcpr))
                    if opts["print"] == True : print(r,"-",c)
                    if opts["print"] == True : print(cpl,"->",nzcpl)
                    if opts["print"] == True : print(cpr,"->",nzcpr)

                    sm1D.at[r,c] = 0.0         
                    for w in range(0,len(allpairs)):
                        li = tuple(allpairs[w][0])
                        ri = tuple(allpairs[w][1])
                        ll = cpl[li]
                        rr = cpr[ri]                
                        #if ll != 0.0 and rr != 0.0 :                    
                        dd = omd[k].at[li,ri]
                        if opts["print"] == True : print("sum : ",ll," | ",rr,"|",dd)
                        #if dd is None : dd = 0.0
                        sm1D.at[r,c] = sm1D.at[r,c] + ll*rr*dd

                    sm1D.at[c,r] = sm1D.at[r,c]
                    if opts["print"] == True : print("\n")

            sm.append(sm1D)
            out = sum(sm)
            out.drop(edge.index[ edge["edge"] == True ],inplace = True,axis=0)
            out.drop(edge.index[ edge["edge"] == True ],inplace = True,axis=1)

        return out
    
    ###
    def _scalar(self):#restituisce la stessa Bspline ma con cp scalari
        sh = shape(dim=self.dim(),codim=1)
        kv = self._kv        
        return Bspline(sh,kv)#cp: default value
    
    ###
    def Galerkin(self,func,opts=None):
        opts = self.prepare_opts(opts)

        sm = self.stiffness_matrix(opts)
        lv = self.load_vector(func,opts)

        sm.replace(np.nan,0.0,inplace=True)
        lv.replace(np.nan,0.0,inplace=True)

        smnp = np.asarray(sm)
        lvnp = np.asarray(lv)

        cp  = np.linalg.solve(smnp,lvnp)
        out = pd.DataFrame(cp,index=sm.index,columns=["cp"])

        for i in range(len(cp)):
            j = out.index[i]
            self._cp[j]  = out.at[j,"cp"]
        return out      
    
    ###
    def approximate(self,func,opts=None):
        #http://hplgit.github.io/INF5620/doc/pub/sphinx-fem/._main_fem002.html
        opts = self.prepare_opts(opts,opts2={"norm":"L1","del-edge":False})

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
                self.set_cp(j, out.iloc[j])
        return out       
    
    ###
    def load_vector(self,func,opts=None): 
        
        opts = self.prepare_opts(opts)

        # escludo subito dal calcolo i termini di bordo 
        #impostando a nan il valore in sm1D
        if opts["del-edge"] == True :
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
        #integrate = lambda *xx : 


        for i in il :

            #controllo che il termine non sia di bordo
            if opts["del-edge"] == True and edge.at[i,"edge"] == True :
                out.at[i,"cp"] = np.nan
                continue

            #i = tuple(i)
            scalar.set_cp(i,1.0)
            ov = self.basis_overlap(i,i,br) #overlap

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
        if opts["del-edge"] == True:
            out = out.drop(edge.index[ edge["edge"] == True ])  
            
        return out       
    
###        
def overlaps(a, b):
    c = [ max(a[0],b[0]) , min(a[1],b[1]) ]
    
    if c[0] >= c[1] :
        return None
    else :
        return c


# In[ ]:




