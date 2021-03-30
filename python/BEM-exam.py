#!/usr/bin/env python
# coding: utf-8

# # Pacchetti

# In[1]:


import pyBspline as Bs
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import scipy
from scipy.misc import derivative
from scipy.optimize import curve_fit
import pandas as pd
#from ipywidgets import interactive
#import ipywidgets as widgets
#from ipywidgets import AppLayout, FloatSlider
#from mpl_toolkits.mplot3d import Axes3D
import copy
import pandas as pd
#from scipy import integrate
#import itertools 
#import time
import os
import scipy.special
import re
import FFT as esFFT
from imp import reload 

###
def norm(x):
    return np.sqrt(np.sum(np.power(x,2.0)))        
      


# In[2]:


#
def plot(fig,n,xB,yB,x,y,c,title,cmap):
    
    ax = fig.add_subplot(n)
    ax.plot(xB, yB, color= "black",label="Bspline")
    sc = ax.scatter(x,y,c=c,cmap=cmap)
    plt.colorbar(sc)
    ax.set_aspect('equal')
    plt.xlim(min(x),max(x))
    plt.ylim(min(y),max(y))
    plt.title(title)
    
    return


# In[35]:


#
def plot_sol_disc(fig,n,xB_arr,yB_arr,x,y,c,title,cmap):
    
    ax = fig.add_subplot(n)
    for xB,yB in zip(xB_arr,yB_arr):
        ax.plot(xB, yB, color= "black",label="Bspline")
    sc = ax.scatter(x,y,c=c,cmap=cmap)
    plt.colorbar(sc)
    ax.set_aspect('equal')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.title(title)
    
    return


# In[3]:


#
def plot_sol(fig,n,xB,yB,x,y,c,title):
    
    ax = fig.add_subplot(n)
    ax.plot(xB, yB, color= "black",label="Bspline")
    sc = ax.scatter(x,y,c=c,cmap=cmap)
    plt.colorbar(sc)
    ax.set_aspect('equal')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.title(title)
    
    return


# In[4]:


#
def plot_matrix(sm,file_png=None):
    #
    sm2 = sm.copy()
    sm2["index"] = sm2.index
    #
    new = sm2.melt(id_vars=['index'])# = sm.index
    #
    new2 = new.copy()
    new2["index"] = [ i[0] for i in new2["index"]]
    new2["variable"] = [ i[0] for i in new2["variable"]]
    #new2["value"] = [np.complex(i) for i in new2["value"] ]
    new2["real"] = np.real(new2["value"])
    new2["imag"] = np.imag(new2["value"])
    new2["abs"] = np.absolute(new2["value"])
    new2["phase"] = np.angle(new2["value"])/np.pi
    df = new2
    
    #
    fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

    cmap = 'RdYlBu'

    #
    ax = fig.add_subplot(221)
    sc = ax.scatter(df["index"],df["variable"],c=df["real"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : real")

    #
    ax = fig.add_subplot(222)
    sc = ax.scatter(df["index"],df["variable"],c=df["imag"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : imag")

    #
    ax = fig.add_subplot(223)
    sc = ax.scatter(df["index"],df["variable"],c=df["abs"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : abs")

    #
    ax = fig.add_subplot(224)
    sc = ax.scatter(df["index"],df["variable"],c=df["phase"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : $\\theta / \pi$")

    plt.tight_layout()
    if file_png is not None :
        plt.savefig(file_png)
    plt.show()
    


# # Basis function

# In[5]:


#definisco la Bspline
sh = Bs.shape(1,2)
#sh.show()

#defiisco i knot vector
P=1 #polinomial degree
N=20 #base caridnality
xminBs = 0.0
xmaxBs = 1.0


#
#kv = Bs.uniform_open_kv(xmin,xmax,p=P,n=N)#Bs.knot_vector(P,N,v)
#kv = periodic_kv(xmin,xmax,p=P,n=N)
kv = Bs.periodic_kv(xminBs,xmaxBs,p=P,n=N)
#kv.show()

#alloco la Bspline
bs = Bs.Bspline(sh,[kv],properties={"periodic":[True]})


# In[5]:


#function
x0 = 0.0
y0 = 0.0
a = 1.0
b = 1.0
def func(t):
    #print(cpz)
    cpx = a*np.cos(2*np.pi*t)+x0#np.random.rand(N)
    cpy = b*np.sin(2*np.pi*t)+y0#np.random.rand(N)
    out = np.zeros(shape=(len(t),2))
    for i in range(len(t)):
        out[i,0] = cpx[i]
        out[i,1] = cpy[i]
    return out


# In[6]:


#control points
t = np.linspace(0,1,N,endpoint=False)
cp = func(t)
for i in range(len(t)):
    #bs._cp[i] = cp[i]
    bs.set_cp(i,cp[i])
cpx = cp[:,0]
cpy = cp[:,1]


# In[7]:


#valutazione della Bspline
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
df = pd.DataFrame(xy)
df = df.rename(columns={0:"x",1:"y"})


# In[8]:


#grafico
fig = plt.figure ( 0 , figsize = ( 10, 8 ) )

bsCopy = bs.copy()
bsCopy.clear_cp()
bsScal = bsCopy._scalar()
br = bs.basis_range()

#
#s = 0.2
ax = fig.add_subplot(111, projection='3d')

for i in range(N):
    print(i+1,"/",N,end="\r")
    u = np.linspace(br.at[(i,),("min",0)],br.at[(i,),("max",0)],100)
    #u = np.linspace(0,1,1000,endpoint=True)
    bsCopy.set_cp(i,bs.get_cp(i))
    bsScal.set_cp(i,1.)
    xyB = bs.evaluate(u)
    zB = bsScal.evaluate(u)
    
    #xyB = xyB[zB != 0]
    #zB  = zB[zB != 0]
    
    #xyB,zB = [ i,j for i,j in zip(xyB,zB) if j != 0.0 ]
    
    ax.plot(xyB[:,0], xyB[:,1],zB,label=str(i))
    
    bsCopy.set_cp(i,[0,0])
    bsScal.set_cp(i,0.)
    
print("Finished")

ax.plot(xy[:,0], xy[:,1],0.0,color="red",label="Bspline")
ax.set_zlim(0,1)
plt.grid(True)
#plt.legend()
plt.show()


# # Triangle

# ## Definition

# In[4]:


#definisco il vettor d'onda
k_in = 20*np.asarray([0.5,np.sqrt(3)/2])
wavevector = np.sqrt(np.sum(np.power(k_in,2.0)))
I = np.complex(0,1)

xmin = -1.5
xmax = 1.5
ymin = -1.5
ymax = 1.5


# In[5]:


#definisco la dimensionaità:
sh = Bs.shape(1,2)
#sh.show()

#defiisco i knot vector
P=1 #polinomial degree
N=100 #base caridnality
xminBs = 0.0
xmaxBs = 1.0


#
#kv = Bs.uniform_open_kv(xmin,xmax,p=P,n=N)#Bs.knot_vector(P,N,v)
#kv = periodic_kv(xmin,xmax,p=P,n=N)
kv = Bs.periodic_kv(xminBs,xmaxBs,p=P,n=N)
#kv.show()

#alloco la Bspline
bs = Bs.Bspline(sh,[kv],properties={"periodic":[True]})


# In[6]:


#files
file_dir = "files/BEM/triangle-periodic/"
suffix = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)+".csv"
suffix_png = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)+".png"


# ## Geometry

# In[7]:


#triangolo
x0 = -0.5
y0 =  -0.5

a = 1.0 / (2+np.sqrt(2))
b = (1.0+np.sqrt(2)) / (2+np.sqrt(2))
delta = b-a

sx = 1.0
sy = 1.0

def triangle_x(i):
    if i <= a :
        return i
    elif i > a and i <= b :
        j = i-a
        return a - j*a/(delta)
    else :
        return 0
    
def triangle_y(i):
    if i <= a :
        return 0
    elif i > a and i <= b :
        j = i-a
        return j*a/(delta)
    else :
        j = i-b
        return a-j
    
def func(t):
    out = np.zeros((len(t),2))
    out[:,0] = [triangle_x(i) for i in t]
    out[:,1] = [triangle_y(i) for i in t]
            
    out[:,0] = sx*out[:,0]/a+x0
    out[:,1] = sy*out[:,1]/a+y0
    return out


# In[8]:


#ATTENZIONE: mi servono dei punti distribuiti in modo uniforme
# per come ho costruito func so che
# func(0.) = func(1.)
# quindi genero un punto in più
t = np.linspace(0,1,N+1,endpoint=True)#[0:-2]
cp = func(t)
for i in range(len(t)):
    #bs._cp[i] = cp[i]
    bs.set_cp(i,cp[i])
cpx = cp[:,0]
cpy = cp[:,1]


# In[9]:


#valutazione della Bspline
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
Txy = func(T)
df = pd.DataFrame(xy)
df = df.rename(columns={0:"x",1:"y"})


# In[10]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#converto in dataframe    
ax = fig.add_subplot(131)
plt.plot(Txy[:,0],Txy[:,1],color="green",label="triangle",linestyle="--")
plt.scatter(cpx,cpy,color="green",label="cp")
#plt.scatter(df["x"], df["y"], color= "red",label="Bspline")
plt.grid()
plt.legend()
ax.set_aspect('equal')

ax = fig.add_subplot(132)
#plt.plot(cpx,cpy,color="green",label="cp",linestyle="--")
#plt.scatter(cpx,cpy,color="green",label="cp")
plt.scatter(df["x"], df["y"], color= "red",label="Bspline")
plt.grid()
plt.legend()
ax.set_aspect('equal')

#real
ax = fig.add_subplot(133)#, projection='3d')
ax.plot(T,df["x"],color="blue",label="x")
ax.plot(T,df["y"],color="green",label="y")
#ax.plot(df["t"],np.real(df["trace"]),color="red",label="trace")
#plt.title("real")
plt.grid(True)
plt.legend()
plt.show()


# ## Stiffness Matrix

# In[11]:


#files
file = file_dir+"stiffness_matrix-n=6-random=False-"+suffix
file


# In[12]:


#stiffness matrix
READ = True
SAVE = False
if os.path.exists(file) and READ == True :
    sm = bs.load("sm-BEM",file)
else :
    sm,out = bs.stiffness_matrix_BEM(k=wavevector,                opts={"print":True,"N":[6],"return_both":True,"ready_sm_BEM":False,"random":False})
    if SAVE == True :
        bs.save("sm-BEM",file)
sm.head()


# In[13]:


#
def plot_matrix(sm,file_png=None):
    #
    sm2 = sm.copy()
    sm2["index"] = sm2.index
    #
    new = sm2.melt(id_vars=['index'])# = sm.index
    #
    new2 = new.copy()
    new2["index"] = [ i[0] for i in new2["index"]]
    new2["variable"] = [ i[0] for i in new2["variable"]]
    #new2["value"] = [np.complex(i) for i in new2["value"] ]
    new2["real"] = np.real(new2["value"])
    new2["imag"] = np.imag(new2["value"])
    new2["abs"] = np.absolute(new2["value"])
    new2["phase"] = np.angle(new2["value"])/np.pi
    df = new2
    
    #
    fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

    cmap = 'RdYlBu'

    #
    ax = fig.add_subplot(221)
    sc = ax.scatter(df["index"],df["variable"],c=df["real"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : real")

    #
    ax = fig.add_subplot(222)
    sc = ax.scatter(df["index"],df["variable"],c=df["imag"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : imag")

    #
    ax = fig.add_subplot(223)
    sc = ax.scatter(df["index"],df["variable"],c=df["abs"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : abs")

    #
    ax = fig.add_subplot(224)
    sc = ax.scatter(df["index"],df["variable"],c=df["phase"],cmap = 'RdYlBu')
    plt.colorbar(sc)
    plt.xlim(min(df["index"]),max(df["index"]))
    plt.ylim(min(df["index"]),max(df["index"]))
    ax.set_aspect('equal')
    plt.title("stiffness matrix : $\\theta / \pi$")

    plt.tight_layout()
    if file_png is not None :
        plt.savefig(file_png)
    plt.show()
    


# In[14]:


#grafico
ile_png = file_dir+"stiffness_matrix-n=6-random=False-"+suffix_png
plot_matrix(sm,file_png)


# ## Single Layer Potential basis

# In[15]:


#files
file = file_dir+"single_layer_potential-"+suffix
file


# In[16]:


#punti XY
Nx= int(xmax-xmin)*20
Ny = int(ymax-ymin)*20
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY0 = np.zeros((Nx*Ny,2))
XY0[:,0] = X.reshape((Nx*Ny,))
XY0[:,1] = Y.reshape((Nx*Ny,))

#tolgo elementi interni
#radius = np.asarray([np.sqrt(np.sum(np.power(i,2.0))) for i in XY])
internal = bs.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XY0[ np.logical_not(internal) ]


# In[17]:


#single layer potentail per le funzioni di base
READ = True
SAVE = True
if os.path.exists(file) and READ == True :
    slp = bs.load("slp-BEM",file)
#else :

# I can update it adding some  points
slp = bs.single_layer_potential_basis_BEM(XY=XY,k=wavevector,                                              opts={"print":True,"N":[6]})
if SAVE == True :
    bs.save("slp-BEM",file)
    
slp.head()


# ## Plane wave

# ### Preparation

# In[18]:


#punti XY
Nx= int(xmax-xmin)*20
Ny = int(ymax-ymin)*20
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY0 = np.zeros((Nx*Ny,2))
XY0[:,0] = X.reshape((Nx*Ny,))
XY0[:,1] = Y.reshape((Nx*Ny,))

#tolgo elementi interni
#radius = np.asarray([np.sqrt(np.sum(np.power(i,2.0))) for i in XY])
internal = bs.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XY0[ np.logical_not(internal) ]


# In[19]:


#plane wave
def plane_wave(xx): # soluzione
    xx = np.asarray(xx)
    theta = np.dot(xx,k_in)
    return np.exp(I*theta)


# In[20]:


#valutazione della Bspline
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
df = pd.DataFrame(xy)
df = df.rename(columns={0:"x",1:"y"})
uinc = plane_wave(xy)


# In[21]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#
#s = 0.2
ax = fig.add_subplot(121, projection='3d')
ax.plot(xy[:,0], xy[:,1],uinc.real,color="blue",label="real")
ax.plot(xy[:,0], xy[:,1],uinc.imag,color="green",label="imag")
ax.plot(xy[:,0], xy[:,1],0.0,color="red",label="Bspline")
plt.grid(True)
plt.legend()

#
ax = fig.add_subplot(122)#, projection='3d')
ax.plot(T,uinc.real,color="blue",label="real")
ax.plot(T,uinc.imag,color="green",label="imag")
plt.grid(True)
plt.legend()

plt.show()


# In[22]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )

Uinc = plane_wave(XY)

cmap = 'RdYlBu'
    
plot(fig,121,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real",cmap)
plot(fig,122,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag",cmap)

plt.show()


# ### Solution

# In[23]:


#files
file_sol = file_dir+"solution-plane_wave-"+suffix
file_lv  = file_dir+"load_vector-plane_wave-"+suffix
file_ind = file_dir+"indirect_solution-plane_wave-"+suffix
print(file_sol)
print(file_lv)
print(file_ind)


# In[24]:


#metodo di Galerkin
READ = False
SAVE = True
if os.path.exists(file_sol) and READ == True :
    sol,Xnp,Valnp = bs.load("sol-BEM",file_sol)
    
if os.path.exists(file_lv) and READ == True :
    lv = bs.load("lv-BEM",file_lv)
    
if os.path.exists(file_ind) and READ == True :
    sol = bs.load("ind_sol-BEM",file_ind)
    
else :
    opts = {"print":True,"ready_sol_BEM":False,"ready_lv_BEM":False,"ready_ind_sol_BEM":False}
    sol,Xnp,Valnp = bs.BEM(uinc=plane_wave,k=wavevector,XY=XY,opts=opts)
    if SAVE == True :
        bs.save("sol-BEM",file_sol)
        bs.save("lv-BEM",file_lv)
        bs.save("ind_sol-BEM",file_ind)
sol.head()


# In[25]:


#uinc
Uinc = plane_wave(XY)#.reshape(Nx,Ny).transpose()
total = Uinc + Valnp


# In[26]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

cmap = 'RdYlBu'

    
plot_sol(fig,331,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real")
plot_sol(fig,334,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag")
plot_sol(fig,337,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Uinc),"$u_{inc}$ : abs")

plot_sol(fig,332,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real,"$u_{scat}$ : real")
plot_sol(fig,335,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag,"$u_{scat}$ : imag")
plot_sol(fig,338,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp),"$u_{scat}$ : abs")

plot_sol(fig,333,df["x"], df["y"],XY[:,0],XY[:,1],total.real,"$u_{tot}$ : real")
plot_sol(fig,336,df["x"], df["y"],XY[:,0],XY[:,1],total.imag,"$u_{tot}$ : imag")
plot_sol(fig,339,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(total),"$u_{tot}$ : abs")

plt.tight_layout()

sol_png = file_dir+"solution-n=6-random=False-"+suffix_png
plt.savefig(sol_png)

plt.show()


# ## Herglotz

# ### Kernel

# In[27]:


#definisco la dimensionaità:
sh = Bs.shape(1,1)
#sh.show()

#defiisco i knot vector
P=0 #polinomial degree
N=12 #base caridnality

#
#kv = Bs.uniform_open_kv(xmin,xmax,p=P,n=N)#Bs.knot_vector(P,N,v)
#kv = periodic_kv(xmin,xmax,p=P,n=N)
kv = Bs.periodic_kv(0.0,2*np.pi,p=P,n=N)
#kv.show()

#alloco la Bspline
kernel = Bs.Bspline(sh,[kv],properties={"periodic":[True],"dtype":np.complex})
#bs.show()

kernel.clear_cp()
kernel.set_cp(1,1)


# In[28]:


#valutazione della Bspline
NN = 1000
T = np.linspace(0.0,2*np.pi,NN,endpoint=False)
y   = kernel.evaluate(T)
#df = pd.DataFrame(xy)
#df = df.rename(columns={0:"x",1:"y"})


# In[29]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#real
ax = fig.add_subplot(111)#, projection='3d')
ax.plot(T,y.real,color="blue",label="real")
ax.plot(T,y.imag,color="green",label="imag")
plt.xlabel(r"$\theta \, \left[ rad \right]$")
plt.title(r"Herglotz kernel $\, g \left( \theta \right) \, : 1 \, $ if $ \,  - \pi/6 < \theta < \pi/6$")
plt.grid(True)
plt.legend()
plt.show()


# In[30]:


#Herglotz
def Herglotz_private(xy,kernel,k,NN):
    xy = np.asarray(xy)
    theta = np.linspace(0.,2*np.pi,NN,endpoint=False)
    cos = np.cos(theta)
    sin = np.sin(theta)
    g = kernel(theta)
    phase = np.outer(xy[:,0],cos) + np.outer(xy[:,1],sin)
    expo = np.exp(1.j*phase*k)
    return np.dot(expo,g)/NN


# In[31]:


#k = 30./4.
NN = 100
def Herglotz(xy): 
    return Herglotz_private(xy,kernel,wavevector,NN)


# In[32]:


#punti XY
Nx= int(xmax-xmin)*20
Ny = int(ymax-ymin)*20
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY = np.zeros((Nx*Ny,2))
XY[:,0] = X.reshape((Nx*Ny,))
XY[:,1] = Y.reshape((Nx*Ny,))

x = XY[:,0]
y = XY[:,1]


# In[33]:


Uinc = Herglotz(XY)


# In[34]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )

cmap = 'RdYlBu'

ax = fig.add_subplot(131)
#ax.plot(xB, yB, color= "black",label="Bspline")
sc = ax.scatter(x,y,c=Uinc.real,cmap=cmap)
plt.colorbar(sc)
ax.set_aspect('equal')
plt.xlim(min(x),max(x))
plt.ylim(min(y),max(y))
plt.title("real")

ax = fig.add_subplot(132)
#ax.plot(xB, yB, color= "black",label="Bspline")
sc = ax.scatter(x,y,c=Uinc.imag,cmap=cmap)
plt.colorbar(sc)
ax.set_aspect('equal')
plt.xlim(min(x),max(x))
plt.ylim(min(y),max(y))
plt.title("imag")

ax = fig.add_subplot(133)
#ax.plot(xB, yB, color= "black",label="Bspline")
sc = ax.scatter(x,y,c=np.absolute(Uinc),cmap=cmap)
plt.colorbar(sc)
ax.set_aspect('equal')
plt.xlim(min(x),max(x))
plt.ylim(min(y),max(y))
plt.title("abs")

#ax = fig.add_subplot(224)
##ax.plot(xB, yB, color= "black",label="Bspline")
#sc = ax.scatter(x,y,c=np.angle(Uinc),cmap=cmap)
#plt.colorbar(sc)
#ax.set_aspect('equal')
#plt.xlim(min(x),max(x))
#plt.ylim(min(y),max(y))
#plt.title("phase")

plt.show()


# ### Preparation

# In[35]:


#punti XY
Nx= int(xmax-xmin)*20
Ny = int(ymax-ymin)*20
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY0 = np.zeros((Nx*Ny,2))
XY0[:,0] = X.reshape((Nx*Ny,))
XY0[:,1] = Y.reshape((Nx*Ny,))

internal = bs.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XY0[ np.logical_not(internal) ]


# In[36]:


#valutazione della Bspline
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
uinc = Herglotz(xy)


# In[37]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#
#s = 0.2
ax = fig.add_subplot(121, projection='3d')
ax.plot(xy[:,0], xy[:,1],uinc.real,color="blue",label="real")
ax.plot(xy[:,0], xy[:,1],uinc.imag,color="green",label="imag")
ax.plot(xy[:,0], xy[:,1],0.0,color="red",label="Bspline")
plt.grid(True)
plt.legend()

#
ax = fig.add_subplot(122)#, projection='3d')
ax.plot(T,uinc.real,color="blue",label="real")
ax.plot(T,uinc.imag,color="green",label="imag")
plt.grid(True)
plt.legend()

#
#ax = fig.add_subplot(133)#, projection='3d')
#ax.plot(T,somma.real,color="blue",label="real")
#ax.plot(T,somma.imag,color="green",label="imag")
#plt.grid(True)
#plt.legend()

plt.show()


# In[38]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )

Uinc = Herglotz(XY)

cmap = 'RdYlBu'
    
plot(fig,121,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real",cmap)
plot(fig,122,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag",cmap)

plt.show()


# ### Solution

# In[39]:


#files
file_sol = file_dir+"solution-Herglotz-"+suffix
file_lv  = file_dir+"load_vector-Herglotz-"+suffix
file_ind = file_dir+"indirect_solution-Herglotz-"+suffix
print(file_sol)
print(file_lv)
print(file_ind)


# In[40]:


#metodo di Galerkin
READ = False
SAVE = True
if os.path.exists(file_sol) and READ == True :
    sol,Xnp,Valnp = bs.load("sol-BEM",file_sol)
    
if os.path.exists(file_lv) and READ == True :
    lv = bs.load("lv-BEM",file_lv)
    
if os.path.exists(file_ind) and READ == True :
    sol = bs.load("ind_sol-BEM",file_ind)
    
else :
    opts = {"print":True,"ready_sol_BEM":False,"ready_lv_BEM":False,"ready_ind_sol_BEM":False}
    sol,Xnp,Valnp = bs.BEM(uinc=Herglotz,k=wavevector,XY=XY,opts=opts)
    if SAVE == True :
        bs.save("sol-BEM",file_sol)
        bs.save("lv-BEM",file_lv)
        bs.save("ind_sol-BEM",file_ind)
sol.head()


# In[41]:


#uinc
Uinc = Herglotz(XY)#.reshape(Nx,Ny).transpose()
total = Uinc + Valnp


# In[42]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

cmap = 'RdYlBu'

    
plot_sol(fig,331,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real")
plot_sol(fig,334,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag")
plot_sol(fig,337,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Uinc),"$u_{inc}$ : abs")

plot_sol(fig,332,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real,"$u_{scat}$ : real")
plot_sol(fig,335,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag,"$u_{scat}$ : imag")
plot_sol(fig,338,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp),"$u_{scat}$ : abs")

plot_sol(fig,333,df["x"], df["y"],XY[:,0],XY[:,1],total.real,"$u_{tot}$ : real")
plot_sol(fig,336,df["x"], df["y"],XY[:,0],XY[:,1],total.imag,"$u_{tot}$ : imag")
plot_sol(fig,339,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(total),"$u_{tot}$ : abs")

plt.tight_layout()

sol_png = file_dir+"solution-n=6-random=False-"+suffix_png
plt.savefig(sol_png)

plt.show()


# # Circle

# ## Definition

# In[27]:


#definissco il vettor d'onda
k_in = 30/4*np.asarray([np.sqrt(3.)/2.,0.5])
wavevector = 7.5#np.sqrt(np.sum(np.power(k_in,2.0)))
print("w:",wavevector)
print("h:",np.pi/(5*wavevector))
I = np.complex(0,1)

xmin = -3
xmax = 3
ymin = -3
ymax = 3


# In[28]:


#definisco la dimensionaità:
sh = Bs.shape(1,2)
#sh.show()

#defiisco i knot vector
P=1 #polinomial degree
N=100 #base caridnality
xminBs = 0.0
xmaxBs = 1.0


#
#kv = Bs.uniform_open_kv(xmin,xmax,p=P,n=N)#Bs.knot_vector(P,N,v)
#kv = periodic_kv(xmin,xmax,p=P,n=N)
kv = Bs.periodic_kv(xminBs,xmaxBs,p=P,n=N)
#kv.show()

#alloco la Bspline
bs = Bs.Bspline(sh,[kv],properties={"periodic":[True]})


# In[29]:


#files
file_dir = "files/BEM/circle/"
suffix = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)+".csv"
suffix_png = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)+".png"


# ## Geometry

# In[18]:


#geometria
x0 = 0.0
y0 = 0.0
a = 1.0
b = 1.0
radius = 1.0
def func(t):
    #print(cpz)
    cpx = a*np.cos(2*np.pi*t)+x0#np.random.rand(N)
    cpy = b*np.sin(2*np.pi*t)+y0#np.random.rand(N)
    out = np.zeros(shape=(len(t),2))
    for i in range(len(t)):
        out[i,0] = cpx[i]
        out[i,1] = cpy[i]
    return out


# In[19]:


#ATTENZIONE: mi servono dei punti distribuiti in modo uniforme
# per come ho costruito func so che
# func(0.) = func(1.)
# quindi genero un punto in più
t = np.linspace(0,1,N+1,endpoint=True)#[0:-2]
cp = func(t)
for i in range(len(t)):
    #bs._cp[i] = cp[i]
    bs.set_cp(i,cp[i])
cpx = cp[:,0]
cpy = cp[:,1]


# In[20]:


#valutazione della Bspline
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
df = pd.DataFrame(xy)
df = df.rename(columns={0:"x",1:"y"})


# In[21]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#converto in dataframe    
ax = fig.add_subplot(131)
#plt.plot(cpx,cpy,color="green",label="cp",linestyle="--")
plt.scatter(cpx,cpy,color="green",label="cp")
#plt.scatter(df["x"], df["y"], color= "red",label="Bspline")
plt.grid()
plt.legend()
ax.set_aspect('equal')

ax = fig.add_subplot(132)
#plt.plot(cpx,cpy,color="green",label="cp",linestyle="--")
#plt.scatter(cpx,cpy,color="green",label="cp")
plt.scatter(df["x"], df["y"], color= "red",label="Bspline")
plt.grid()
plt.legend()
ax.set_aspect('equal')

#real
ax = fig.add_subplot(133)#, projection='3d')
ax.plot(T,df["x"],color="blue",label="x")
ax.plot(T,df["y"],color="green",label="y")
#ax.plot(df["t"],np.real(df["trace"]),color="red",label="trace")
#plt.title("real")
plt.grid(True)
plt.legend()
plt.show()


# In[22]:


filename = file_dir+"control_points-"+suffix
a = bs.save("cp",filename)
#bs.load("cp",filename)


# ## Stiffness Matrix

# In[10]:


#stiffness matrix
READ = True
SAVE = False
file = file_dir+"stiffness_matrix-n=6-"+suffix
print(file)

if os.path.exists(file) and READ == True :
    sm = bs.load("sm-BEM",file)
else :
    sm,out = bs.stiffness_matrix_BEM(k=wavevector,                                 opts={"print":True,"N":[6],"ready_sm_BEM":False,"return_both":True})
    if SAVE == True :
        bs.save("sm-BEM",file)
sm.head()


# In[17]:


#grafico
file_png = file_dir+"stiffness_matrix-n=6-random=False-"+suffix_png
plot_matrix(sm,file_png)


# ## Single Layer Potential basis

# In[18]:


#punti XY
Nx= int(xmax-xmin)*10
Ny = int(ymax-ymin)*10
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY0 = np.zeros((Nx*Ny,2))
XY0[:,0] = X.reshape((Nx*Ny,))
XY0[:,1] = Y.reshape((Nx*Ny,))

internal = bs.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XY0[ np.logical_not(internal) ]
print(len(XY))#," = ",len(XYslp)/3600,"h")

#tolgo elementi interni
#radius = np.asarray([np.sqrt(np.sum(np.power(i,2.0))) for i in XY])
#XYslp = XY#[radius > 1.0]
#print(len(XYslp)," = ",len(XYslp)*3/3600,"h")


# In[19]:


#single layer potential per le funzioni di base
READ = True
SAVE = True

file = file_dir+"single_layer_potential-"+suffix
print(file)

if os.path.exists(file) and READ == True :
    slp = bs.load("slp-BEM",file)
#else :

# I can update it adding some  points
slp = bs.single_layer_potential_basis_BEM(XY=XY,k=wavevector,                                              opts={"print":True,"N":[6]})
if SAVE == True :
    bs.save("slp-BEM",file)
    
slp.head()


# ## Plane wave

# ### Preparation

# In[20]:


#punti XY
Nx= int(xmax-xmin)*10
Ny = int(ymax-ymin)*10
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY0 = np.zeros((Nx*Ny,2))
XY0[:,0] = X.reshape((Nx*Ny,))
XY0[:,1] = Y.reshape((Nx*Ny,))

internal = bs.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XY0[ np.logical_not(internal) ]


# In[21]:


#plane_wave
def plane_wave(xx): # soluzione
    xx = np.asarray(xx)
    theta = np.dot(xx,k_in)
    return np.exp(I*theta)


# In[22]:


#valutazione della Bspline
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
uinc = plane_wave(xy)


# In[23]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#
#s = 0.2
ax = fig.add_subplot(121, projection='3d')
ax.plot(xy[:,0], xy[:,1],uinc.real,color="blue",label="real")
ax.plot(xy[:,0], xy[:,1],uinc.imag,color="green",label="imag")
ax.plot(xy[:,0], xy[:,1],0.0,color="red",label="Bspline")
plt.grid(True)
plt.legend()

#
ax = fig.add_subplot(122)#, projection='3d')
ax.plot(T,uinc.real,color="blue",label="real")
ax.plot(T,uinc.imag,color="green",label="imag")
plt.grid(True)
plt.legend()

#
#ax = fig.add_subplot(133)#, projection='3d')
#ax.plot(T,somma.real,color="blue",label="real")
#ax.plot(T,somma.imag,color="green",label="imag")
#plt.grid(True)
#plt.legend()

plt.show()


# In[24]:


#uinc
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )

Uinc = plane_wave(XY)

cmap = 'RdYlBu'
    
plot(fig,121,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real",cmap)
plot(fig,122,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag",cmap)

plt.show()


# ### Solution

# In[82]:


#metodo di Galerkin
READ = True
SAVE = False

#
file_sol = file_dir+"solution-plane_wave-"+suffix
file_lv  = file_dir+"load_vector-plane_wave-"+suffix
file_ind = file_dir+"indirect_solution-plane_wave-"+suffix
print(file_sol)
print(file_lv)
print(file_ind)

if os.path.exists(file_sol) and READ == True :
    sol,Xnp,Valnp = bs.load("sol-BEM",file_sol)
    
if os.path.exists(file_lv) and READ == True :
    lv = bs.load("lv-BEM",file_lv)
    
if os.path.exists(file_ind) and READ == True :
    sol = bs.load("ind_sol-BEM",file_ind)
    
else :
    opts = {"print":True,"ready_sol_BEM":False,"ready_lv_BEM":False,"ready_ind_sol_BEM":False}
    sol,Xnp,Valnp = bs.BEM(uinc=plane_wave,k=wavevector,XY=XY,opts=opts)
    if SAVE == True :
        bs.save("sol-BEM",file_sol)
        bs.save("lv-BEM",file_lv)
        bs.save("ind_sol-BEM",file_ind)
sol.head()


# In[27]:


#uinc
Uinc = plane_wave(XY)#.reshape(Nx,Ny).transpose()
total = Uinc + Valnp


# In[28]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

cmap = 'RdYlBu'

   
plot_sol(fig,331,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real")
plot_sol(fig,334,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag")
plot_sol(fig,337,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Uinc),"$u_{inc}$ : abs")

plot_sol(fig,332,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real,"$u_{scat}$ : real")
plot_sol(fig,335,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag,"$u_{scat}$ : imag")
plot_sol(fig,338,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp),"$u_{scat}$ : abs")

plot_sol(fig,333,df["x"], df["y"],XY[:,0],XY[:,1],total.real,"$u_{tot}$ : real")
plot_sol(fig,336,df["x"], df["y"],XY[:,0],XY[:,1],total.imag,"$u_{tot}$ : imag")
plot_sol(fig,339,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(total),"$u_{tot}$ : abs")

plt.tight_layout()

sol_png = file_dir+"solution-n=6-random=False-"+suffix_png
plt.savefig(sol_png)

plt.show()


# ### Analytic solution

# In[78]:


#
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
uinc = plane_wave(xy)


# In[79]:


#
fft = esFFT.FFT(uinc,opts={"plot":True})

#
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )
ax = fig.add_subplot(111)
plt.plot(fft.index, np.real(fft["fft"]),color="blue" ,label="real")#,marker="+")
plt.plot(fft.index, np.imag(fft["fft"]),color="green",label="imag")#,marker="x")
plt.xlim(-20,20)
plt.legend()
plt.grid(True)
plt.title("Fourier Transform")
plt.show()
#
fft.head()


# In[80]:


#
out,analytic = esFFT.analytic_solution_circle(uinc,XY,wmin=-15,wmax=15,radius=radius,                                              wavevector=wavevector,opts={"return":"both"})

analytic_tot = Uinc + analytic

out.head()


# In[83]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

cmap = 'RdYlBu'
    
plot_sol(fig,331,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real,"$u_{scat}^{Galerkin}$ : real")
plot_sol(fig,334,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag,"$u_{scat}^{Galerkin}$ : imag")
plot_sol(fig,337,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp),"$u_{scat}^{Galerkin}$ : abs")

plot_sol(fig,332,df["x"], df["y"],XY[:,0],XY[:,1],analytic.real,"$u_{scat}^{analytic}$ : real")
plot_sol(fig,335,df["x"], df["y"],XY[:,0],XY[:,1],analytic.imag,"$u_{scat}^{analytic}$ : imag")
plot_sol(fig,338,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(analytic),"$u_{scat}^{analytic}$ : abs")

plot_sol(fig,333,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real-analytic.real,         "$u_{scat}^{Galerkin}-u_{scat}^{analytic}$ : real")
plot_sol(fig,336,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag-analytic.imag,         "$u_{scat}^{Galerkin}-u_{scat}^{analytic}$ : imag")
plot_sol(fig,339,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp-analytic),         "$u_{scat}^{Galerkin}-u_{scat}^{analytic}$ : abs")


plt.tight_layout()

sol_png = file_dir+"solution-analytic-n=6-random=False-"+suffix_png
plt.savefig(sol_png)

plt.show()


# ## Herglotz

# ### Kernel

# In[29]:


#definisco la dimensionaità:
sh = Bs.shape(1,1)
#sh.show()

#defiisco i knot vector
P=0 #polinomial degree
N=12 #base caridnality

#
#kv = Bs.uniform_open_kv(xmin,xmax,p=P,n=N)#Bs.knot_vector(P,N,v)
#kv = periodic_kv(xmin,xmax,p=P,n=N)
kv = Bs.periodic_kv(0.0,2*np.pi,p=P,n=N)
#kv.show()

#alloco la Bspline
kernel = Bs.Bspline(sh,[kv],properties={"periodic":[True],"dtype":np.complex})
#bs.show()

kernel.clear_cp()
kernel.set_cp(0,1)


# In[30]:


#
NN = 1000
T = np.linspace(0.0,2*np.pi,NN,endpoint=False)
y   = kernel.evaluate(T)
#df = pd.DataFrame(xy)
#df = df.rename(columns={0:"x",1:"y"})


# In[31]:


#
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#real
ax = fig.add_subplot(111)#, projection='3d')
ax.plot(T,y.real,color="blue",label="real")
ax.plot(T,y.imag,color="green",label="imag")
plt.xlabel(r"$\theta \, \left[ rad \right]$")
plt.title(r"Herglotz kernel $\, g \left( \theta \right) \, : 1 \, $ if $ \,  - \pi/6 < \theta < \pi/6$")
plt.grid(True)
plt.legend()
plt.show()


# In[32]:


#
def Herglotz_private(xy,kernel,k,NN):
    xy = np.asarray(xy)
    theta = np.linspace(0.,2*np.pi,NN,endpoint=False)
    cos = np.cos(theta)
    sin = np.sin(theta)
    g = kernel(theta)
    phase = np.outer(xy[:,0],cos) + np.outer(xy[:,1],sin)
    expo = np.exp(1.j*phase*k)
    return np.dot(expo,g)/NN


# In[33]:


#k = 30./4.
NN = 100
def Herglotz(xy): 
    return Herglotz_private(xy,kernel,wavevector,NN)


# In[34]:


#
Nx= int(xmax-xmin)*10
Ny = int(ymax-ymin)*10
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY = np.zeros((Nx*Ny,2))
XY[:,0] = X.reshape((Nx*Ny,))
XY[:,1] = Y.reshape((Nx*Ny,))

x = XY[:,0]
y = XY[:,1]


# In[35]:


Uinc = Herglotz(XY)


# In[36]:


#
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )

cmap = 'RdYlBu'

ax = fig.add_subplot(131)
#ax.plot(xB, yB, color= "black",label="Bspline")
sc = ax.scatter(x,y,c=Uinc.real,cmap=cmap)
plt.colorbar(sc)
ax.set_aspect('equal')
plt.xlim(min(x),max(x))
plt.ylim(min(y),max(y))
plt.title("real")

ax = fig.add_subplot(132)
#ax.plot(xB, yB, color= "black",label="Bspline")
sc = ax.scatter(x,y,c=Uinc.imag,cmap=cmap)
plt.colorbar(sc)
ax.set_aspect('equal')
plt.xlim(min(x),max(x))
plt.ylim(min(y),max(y))
plt.title("imag")

ax = fig.add_subplot(133)
#ax.plot(xB, yB, color= "black",label="Bspline")
sc = ax.scatter(x,y,c=np.absolute(Uinc),cmap=cmap)
plt.colorbar(sc)
ax.set_aspect('equal')
plt.xlim(min(x),max(x))
plt.ylim(min(y),max(y))
plt.title("abs")

#ax = fig.add_subplot(224)
##ax.plot(xB, yB, color= "black",label="Bspline")
#sc = ax.scatter(x,y,c=np.angle(Uinc),cmap=cmap)
#plt.colorbar(sc)
#ax.set_aspect('equal')
#plt.xlim(min(x),max(x))
#plt.ylim(min(y),max(y))
#plt.title("phase")

plt.show()


# ### Preparation

# In[95]:


#
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
uinc = Herglotz(xy)


# In[37]:


#
Nx= int(xmax-xmin)*10
Ny = int(ymax-ymin)*10
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY = np.zeros((Nx*Ny,2))
XY[:,0] = X.reshape((Nx*Ny,))
XY[:,1] = Y.reshape((Nx*Ny,))

internal = bs.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XY0[ np.logical_not(internal) ]

x = XY[:,0]
y = XY[:,1]


# In[38]:


#
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )

#
#s = 0.2
ax = fig.add_subplot(121, projection='3d')
ax.plot(xy[:,0], xy[:,1],uinc.real,color="blue",label="real")
ax.plot(xy[:,0], xy[:,1],uinc.imag,color="green",label="imag")
ax.plot(xy[:,0], xy[:,1],0.0,color="red",label="Bspline")
plt.grid(True)
plt.legend()

#
ax = fig.add_subplot(122)#, projection='3d')
ax.plot(T,uinc.real,color="blue",label="real")
ax.plot(T,uinc.imag,color="green",label="imag")
plt.grid(True)
plt.legend()

#
#ax = fig.add_subplot(133)#, projection='3d')
#ax.plot(T,somma.real,color="blue",label="real")
#ax.plot(T,somma.imag,color="green",label="imag")
#plt.grid(True)
#plt.legend()

plt.show()


# In[39]:


#
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )

Uinc = Herglotz(XY)

cmap = 'RdYlBu'
    
plot(fig,121,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real",cmap)
plot(fig,122,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag",cmap)

plt.show()


# ### Solution

# In[40]:


#
READ = True
SAVE = False

file_sol = file_dir+"solution-Herglotz-"+suffix
file_lv  = file_dir+"load_vector-Herglotz-"+suffix
file_ind = file_dir+"indirect_solution-Herglotz-"+suffix
print(file_sol)
print(file_lv)
print(file_ind)

if os.path.exists(file_sol) and READ == True :
    sol,Xnp,Valnp = bs.load("sol-BEM",file_sol)
    
if os.path.exists(file_lv) and READ == True :
    lv = bs.load("lv-BEM",file_lv)
    
if os.path.exists(file_ind) and READ == True :
    sol = bs.load("ind_sol-BEM",file_ind)
    
else :
    opts = {"print":True,"ready_sol_BEM":False,"ready_lv_BEM":False,"ready_ind_sol_BEM":False}
    sol,Xnp,Valnp = bs.BEM(uinc=Herglotz,k=wavevector,XY=XY,opts=opts)
    if SAVE == True :
        bs.save("sol-BEM",file_sol)
        bs.save("lv-BEM",file_lv)
        bs.save("ind_sol-BEM",file_ind)
sol.head()


# In[41]:


#
Uinc = Herglotz(XY)#.reshape(Nx,Ny).transpose()
total = Uinc + Valnp


# In[42]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

cmap = 'RdYlBu'
    
plot_sol(fig,331,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real")
plot_sol(fig,334,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag")
plot_sol(fig,337,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Uinc),"$u_{inc}$ : abs")

plot_sol(fig,332,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real,"$u_{scat}$ : real")
plot_sol(fig,335,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag,"$u_{scat}$ : imag")
plot_sol(fig,338,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp),"$u_{scat}$ : abs")

plot_sol(fig,333,df["x"], df["y"],XY[:,0],XY[:,1],total.real,"$u_{tot}$ : real")
plot_sol(fig,336,df["x"], df["y"],XY[:,0],XY[:,1],total.imag,"$u_{tot}$ : imag")
plot_sol(fig,339,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(total),"$u_{tot}$ : abs")

plt.tight_layout()

sol_png = file_dir+"solution-n=6-random=False-"+suffix_png
plt.savefig(sol_png)

plt.show()


# ### Analytic solution

# In[43]:


#
NN = 1000
T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
xy   = bs.evaluate(T)
uinc = Herglotz(xy)


# In[44]:


#
fft = esFFT.FFT(uinc,opts={"plot":True})

#
fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )
ax = fig.add_subplot(111)
plt.plot(fft.index, np.real(fft["fft"]),color="blue" ,label="real")#,marker="+")
plt.plot(fft.index, np.imag(fft["fft"]),color="green",label="imag")#,marker="x")
plt.xlim(-20,20)
plt.legend()
plt.grid(True)
plt.title("Fourier Transform")
plt.show()
#
fft.head()


# In[75]:


#
out,analytic = esFFT.analytic_solution_circle(uinc,XY,wmin=-15,wmax=15,radius=radius,                                              wavevector=wavevector,opts={"return":"both"})

analytic_tot = Uinc + analytic

out.head()


# In[77]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )

cmap = 'RdYlBu'
    
plot_sol(fig,331,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real,"$u_{scat}^{Galerkin}$ : real")
plot_sol(fig,334,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag,"$u_{scat}^{Galerkin}$ : imag")
plot_sol(fig,337,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp),"$u_{scat}^{Galerkin}$ : abs")

plot_sol(fig,332,df["x"], df["y"],XY[:,0],XY[:,1],analytic.real,"$u_{scat}^{analytic}$ : real")
plot_sol(fig,335,df["x"], df["y"],XY[:,0],XY[:,1],analytic.imag,"$u_{scat}^{analytic}$ : imag")
plot_sol(fig,338,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(analytic),"$u_{scat}^{analytic}$ : abs")

plot_sol(fig,333,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real-analytic.real,         "$u_{scat}^{Galerkin}-u_{scat}^{analytic}$ : real")
plot_sol(fig,336,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag-analytic.imag,         "$u_{scat}^{Galerkin}-u_{scat}^{analytic}$ : imag")
plot_sol(fig,339,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp-analytic),         "$u_{scat}^{Galerkin}-u_{scat}^{analytic}$ : abs")


plt.tight_layout()

sol_png = file_dir+"solution-analytic-n=6-random=False-"+suffix_png
plt.savefig(sol_png)

plt.show()


# # Two circles

# ## Definition

# In[5]:


#definisco la dimensionaità:
sh = Bs.shape(1,2)
#sh.show()

#defiisco i knot vector
P=1 #polinomial degree
N=20 #base caridnality
xminBs = 0.0
xmaxBs = 1.0


#
#kv = Bs.uniform_open_kv(xmin,xmax,p=P,n=N)#Bs.knot_vector(P,N,v)
#kv = periodic_kv(xmin,xmax,p=P,n=N)
kv = Bs.periodic_kv(xminBs,xmaxBs,p=P,n=N)
#kv.show()

#alloco la Bspline
bs = Bs.Bspline(sh,[kv],properties={"periodic":[True]})


# In[6]:


#
k_in = 30/4*np.asarray([np.sqrt(3.)/2.,0.5])
wavevector = 7.5#np.sqrt(np.sum(np.power(k_in,2.0)))
print("w:",wavevector)
print("h:",np.pi/(5*wavevector))
I = np.complex(0,1)

xmin = -3
xmax = 3
ymin = -3
ymax = 3


# In[7]:


#function
x0 = 0.0
y0 = 0.0
a = 1.0
b = 1.0
radius = 1.0
def func(t):
    #print(cpz)
    cpx = a*np.cos(2*np.pi*t)+x0#np.random.rand(N)
    cpy = b*np.sin(2*np.pi*t)+y0#np.random.rand(N)
    out = np.zeros(shape=(len(t),2))
    for i in range(len(t)):
        out[i,0] = cpx[i]
        out[i,1] = cpy[i]
    return out


# In[8]:


#ATTENZIONE: mi servono dei punti distribuiti in modo uniforme
# per come ho costruito func so che
# func(0.) = func(1.)
# quindi genero un punto in più
t = np.linspace(0,1,N+1,endpoint=True)#[0:-2]
cp = func(t)
for i in range(len(t)):
    #bs._cp[i] = cp[i]
    bs.set_cp(i,cp[i])
cpx = cp[:,0]
cpy = cp[:,1]


# In[9]:


#bs.control_points()


# In[10]:


#files
file_dir = "files/BEM/circle/"
suffix = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)#+".csv"
suffix_png = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)+".png"

directory = "files/BEM/two-circles/"


# ## Stiffness Matrix

# In[11]:


#stiffness matrix
READ = True
SAVE = False
file = directory+"bs-sm-BEM-"+suffix+".csv"
print(file)

if os.path.exists(file) and READ == True :
    sm = bs.load("sm-BEM",file)
else :
    sm = bs.stiffness_matrix_BEM(k=wavevector,opts={"print":True,"N":[6],"ready_sm_BEM":False})
    if SAVE == True :
        bs.save("sm-BEM",file)
sm.head()


# In[12]:


#grafico
sm = bs.stiffness_matrix_BEM()
plot_matrix(sm)


# ## Single Layer Potential Basis

# In[39]:


#punti XY
Nx= 80
Ny = 80
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY0 = np.zeros((Nx*Ny,2))
XY0[:,0] = X.reshape((Nx*Ny,))
XY0[:,1] = Y.reshape((Nx*Ny,))

internal = bs.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XY0[ np.logical_not(internal) ]
print(len(XY))#," = ",len(XYslp)/3600,"h")

#tolgo elementi interni
#radius = np.asarray([np.sqrt(np.sum(np.power(i,2.0))) for i in XY])
#XYslp = XY#[radius > 1.0]
#print(len(XYslp)," = ",len(XYslp)*3/3600,"h")


# In[41]:


#Single Layer Potential per le funzioni di base
READ = True
SAVE = True

file = directory+"bs-slpB-BEM-"+suffix+".csv"
print(file)

if os.path.exists(file) and READ == True :
    slp = bs.load("slp-BEM",file)
#else :

# I can update it adding some  points
slpB = bs.single_layer_potential_basis_BEM(XY=XY,k=wavevector,opts={"print":True,"N":[6]})

if SAVE == True :
    bs.save("slp-BEM",file)
    
slp.head()


# In[42]:


#traslo
deltaX = int(1./(x[1]-x[0]))*(x[1]-x[0])
deltaY = int(1./(y[1]-y[0]))*(y[1]-y[0])

rL = [-deltaX,deltaY]
rR = [ deltaX, -deltaY]


# In[43]:


#valutazione della Bspline
NN = 1000
t = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
df = pd.DataFrame(index=np.arange(0,len(t)),columns=np.arange(0,6),dtype=object)
index = [  ("left","x") , ("left","y"),          ("right","x") , ("right","y")]
mi = pd.MultiIndex.from_tuples(index)
df = df.reindex(columns=mi)

df["t"]      = t

xy = bs.evaluate(t)
#c = center.evaluate(t)
#d =  down.evaluate(t)

df[("left","x")]  = xy[:,0] + rL[0]
df[("left","y")]  = xy[:,1] + rL[1]

df[("right","x")] = xy[:,0] + rR[0]
df[("right","y")] = xy[:,1] + rR[1]


df


# In[44]:


# due cerchi
Left = bs.copy()
#print(Left._slp_BEM.shape)
Left.traslate_cp(rL)
#print(Left._slp_BEM.shape)
#print(Left._ready_slp_BEM)
Left.traslate_slpB(rL)
#print(Left._slp_BEM.shape)
#print(Left._ready_slp_BEM)

Right = bs.copy()
Right.traslate_cp(rR)
Right.traslate_slpB(rR)


# In[46]:


# punti interni
XYleft = np.asarray([ np.asarray(i) for i in Left._slp_BEM.index ])
XYright = np.asarray([ np.asarray(i) for i in Right._slp_BEM.index ])

internal = Left.internal_points(XY=XY0,NN=1000,xmin=0.,xmax=1.,opts=None)
XYL = XY0[ np.logical_not(internal) ]

internal = Right.internal_points(XY=XYL,NN=1000,xmin=0.,xmax=1.,opts=None)
XY = XYL[ np.logical_not(internal)] 
print(len(XY))


# In[47]:


#grafico
fig = plt.figure ( 0 , figsize = ( 15, 5 ) )
#
ax = fig.add_subplot(111)
plt.plot(df[("left","x")], df[("left","y")], color= "black",label="up")
plt.plot(df[("right","x")], df[("right","y")], color= "black",label="center")
plt.scatter(XY[:,0],XY[:,1], color= "blue",label="points",s=0.1)

#plt.scatter(XYleft[:,0],XYleft[:,1], color= "red",label="points",s=0.1)
#plt.scatter(XYright[:,0],XYright[:,1], color= "green",label="points",s=0.1)

#plt.grid()
#plt.legend()
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
ax.set_aspect('equal')
plt.show()


# In[48]:


#Single Layer Potential per le funzioni di base
READ = True
SAVE = True

file = directory+"Left-slpB-BEM-"+suffix+".csv"
print(file)

if os.path.exists(file) and READ == True :
    slpBL = Left.load("slp-BEM",file)
#else :

slpBL = Left.single_layer_potential_basis_BEM(XY=XY,k=wavevector,                                              opts={"print":True,"prec":1e-6,"ready_slp_BEM":True})
if SAVE == True :
    Left.save("slp-BEM",directory+"Left-slpB-BEM-"+suffix+".csv")
    
slpBL.head()


# In[49]:


#Single Layer Potential per le funzioni di base
READ = True
SAVE = True

file = directory+"Right-slpB-BEM-"+suffix+".csv"
print(file)

if os.path.exists(file) and READ == True :
    slpBR = Right.load("slp-BEM",file)
#else :

slpBR = Right.single_layer_potential_basis_BEM(XY=XY,k=wavevector,                                              opts={"print":True,"prec":1e-6,"ready_slp_BEM":True})
if SAVE == True :
    Right.save("slp-BEM",directory+"Right-slpB-BEM-"+suffix+".csv")
    
slpBR.head()


# ## Solution

# In[50]:


#plane wave
def plane_wave(xx): # soluzione
    xx = np.asarray(xx)
    theta = np.dot(xx,k_in)
    return np.exp(1.j*theta)


# In[51]:


# stiffness matrix
SM = pd.read_csv(directory+"SM-"+suffix+".csv",skiprows=2)
del SM["index"]
SM = SM.applymap(np.complex)
index,mi,last,first = Bs.multi_index_disconnected([Left,Right],[0,1],None)
df0 = pd.DataFrame(data = np.asarray(SM),columns=index,index=index)
SM = df0.reindex(columns=mi,index=mi)
Bs._sm_BEM_disc = SM.copy()
Bs._ready_sm_BEM_disc = True
SM.head()


# In[ ]:


#metodo di Galerkin
bs_arr = [Left,Right]
opts = {"print":True,"N":[6]}
SLP,Xnp,Valnp = Bs.BEM_disconnected(bs=bs_arr,uinc=plane_wave,k=wavevector,XY=XY,opts=opts)
SLP.head()


# In[53]:


SLP.to_csv(directory+"solution-"+suffix+".csv",index_label="index")
Bs._sm_BEM_disc.to_csv(directory+"SM-"+suffix+".csv",index_label="index")
Bs._lv_BEM_disc.to_csv(directory+"LV-"+suffix+".csv",index_label="index")


# In[ ]:


plot_matrix(Bs.get_block(SM,0,0,drop=True))


# In[ ]:


plot_matrix(Bs.get_block(SM,1,1,drop=True))


# In[ ]:


plot_matrix(Bs.get_block(SM,0,1,drop=True))


# In[ ]:


plot_matrix(Bs.get_block(SM,1,0,drop=True))


# In[54]:


#uinc
Uinc = plane_wave(XY)#.reshape(Nx,Ny).transpose()
scat = np.asarray(SLP["value"],dtype=np.complex)#Valnp
total = Uinc + scat


# In[55]:


#grafico
fig = plt.figure ( 0 , figsize = ( 13 , 13 ) )

cmap = 'RdYlBu'

x = [df[("left","x")],df[("right","x")]]
y = [df[("left","y")],df[("right","y")]]
   
plot_sol_disc(fig,331,x,y,XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real",cmap)
plot_sol_disc(fig,334,x,y,XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag",cmap)
plot_sol_disc(fig,337,x,y,XY[:,0],XY[:,1],np.absolute(Uinc),"$u_{inc}$ : abs",cmap)

plot_sol_disc(fig,332,x,y,XY[:,0],XY[:,1],scat.real,"$u_{scat}$ : real",cmap)
plot_sol_disc(fig,335,x,y,XY[:,0],XY[:,1],scat.imag,"$u_{scat}$ : imag",cmap)
plot_sol_disc(fig,338,x,y,XY[:,0],XY[:,1],np.absolute(scat),"$u_{scat}$ : abs",cmap)

plot_sol_disc(fig,333,x,y,XY[:,0],XY[:,1],total.real,"$u_{tot}$ : real",cmap)
plot_sol_disc(fig,336,x,y,XY[:,0],XY[:,1],total.imag,"$u_{tot}$ : imag",cmap)
plot_sol_disc(fig,339,x,y,XY[:,0],XY[:,1],np.absolute(total),"$u_{tot}$ : abs",cmap)

plt.tight_layout()

sol_png = directory+"solution-2-"+suffix_png
plt.savefig(sol_png)

plt.show()


# # Convergence

# In[5]:


#vettor d'onda
k_in = 30/4*np.asarray([np.sqrt(3.)/2.,0.5])
wavevector = 7.5#np.sqrt(np.sum(np.power(k_in,2.0)))
print("w:",wavevector)
print("h:",np.pi/(5*wavevector))
I = np.complex(0,1)

xmin = -3
xmax = 3
ymin = -3
ymax = 3


# ## Preparation

# In[6]:


#function
x0 = 0.0
y0 = 0.0
a = 1.0
b = 1.0
radius = 1.0
def func(t):
    #print(cpz)
    cpx = a*np.cos(2*np.pi*t)+x0#np.random.rand(N)
    cpy = b*np.sin(2*np.pi*t)+y0#np.random.rand(N)
    out = np.zeros(shape=(len(t),2))
    for i in range(len(t)):
        out[i,0] = cpx[i]
        out[i,1] = cpy[i]
    return out


# In[7]:


#punti XY
Nx= 40
Ny = 40
x = np.linspace(xmin,xmax,Nx)
y = np.linspace(ymin,ymax,Ny)
X,Y = np.meshgrid(x,y)

XY0 = np.zeros((Nx*Ny,2))
XY0[:,0] = X.reshape((Nx*Ny,))
XY0[:,1] = Y.reshape((Nx*Ny,))
print(len(XY0))


# In[8]:


#plane wave
def plane_wave(xx): # soluzione
    xx = np.asarray(xx)
    theta = np.dot(xx,k_in)
    return np.exp(I*theta)


# In[13]:


#definisco la dimensionaità:
sh = Bs.shape(1,2)
#sh.show()

#defiisco i knot vector
P=1 #polinomial degree
    
xminBs = 0.0
xmaxBs = 1.0
    
#N_arr = [10,20,30,40,50,60,70,80,90,100]
#N_arr = np.arange(11,20)
N_arr = list(np.arange(10,21))+list(np.arange(30,110,10))
#N_arr = [10]


# ## Cycle

# In[ ]:


#ciclo
j = 0
n = 4
for N in N_arr:    
    
    count = str(j)+"/"+str(len(N_arr))
    
    print(count," : preparazione")
    
    #
    file_dir = "files/BEM/circle-convergence/n="+str(n)+"/"
    suffix = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)+"-n="+str(n)
    
    
    kv = Bs.periodic_kv(xminBs,xmaxBs,p=P,n=N)
    #kv.show()

    #alloco la Bspline
    bs = Bs.Bspline(sh,[kv],properties={"periodic":[True]})
    
    t = np.linspace(0,1,N+1,endpoint=True)#[0:-2]
    cp = func(t)
    for i in range(len(t)):
        bs.set_cp(i,cp[i])
        
    #
    internal = bs.internal_points(XY=XY0,NN=1000,xmin=xminBs,xmax=xmaxBs,opts=None)
    XY = XY0[ np.logical_not(internal) ]
    
    #bs.load("sm-BEM"     ,file_dir+"stiffness_matrix-"            +suffix+".csv")    
    #bs.load("slp-BEM"    ,file_dir+"single_layer_potential-"      +suffix+".csv")    
    #bs.load("sol-BEM"    ,file_dir+"solution-plane_wave-"         +suffix+".csv")
    #bs.load("ind_sol-BEM",file_dir+"indirect_solution-plane_wave-"+suffix+".csv")
    
    #
    print(count," : BEM")
    opts={"print":True,"N":[n]} 
    sol,Xnp,Valnp = bs.BEM(uinc=plane_wave,k=wavevector,XY=XY,opts=opts)

    print(count," : saving")
    bs.save("sm-BEM"     ,file_dir+"stiffness_matrix-"            +suffix+".csv")    
    bs.save("slp-BEM"    ,file_dir+"single_layer_potential-"      +suffix+".csv")    
    bs.save("sol-BEM"    ,file_dir+"solution-plane_wave-"         +suffix+".csv")
    #bs.save("lv-BEM"     ,file_dir+"load_vector-plane_wave-"      +suffix+".csv")
    bs.save("ind_sol-BEM",file_dir+"indirect_solution-plane_wave-"+suffix+".csv")
        
    Uinc = plane_wave(XY)
    total = Uinc + Valnp
    
    index = [tuple(i) for i in XY]
    grafico = pd.DataFrame(index=index,columns=["xy","inc","scat","tot"],dtype=object)
    grafico["xy"] = [ i for i in XY]
    grafico["inc"] = Uinc 
    grafico["scat"] = Valnp
    grafico["tot"] = total

    grafico.to_csv(file_dir+"grafico-plane_wave-"+suffix+".csv",index_label="index")
    
    #grafico
    print(count," : grafico")
    fig = plt.figure ( 0 , figsize = ( 15 , 10 ) )  
    
    cmap = 'RdYlBu'
    
    NN = 1000
    T = np.linspace(xminBs,xmaxBs,NN,endpoint=True)
    xy   = bs.evaluate(T)
    df = pd.DataFrame(xy)
    df = df.rename(columns={0:"x",1:"y"})
   
    plot_sol(fig,331,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.real,"$u_{inc}$ : real",cmap)
    plot_sol(fig,334,df["x"], df["y"],XY[:,0],XY[:,1],Uinc.imag,"$u_{inc}$ : imag",cmap)
    plot_sol(fig,337,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Uinc),"$u_{inc}$ : abs",cmap)

    plot_sol(fig,332,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.real,"$u_{scat}$ : real",cmap)
    plot_sol(fig,335,df["x"], df["y"],XY[:,0],XY[:,1],Valnp.imag,"$u_{scat}$ : imag",cmap)
    plot_sol(fig,338,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(Valnp),"$u_{scat}$ : abs",cmap)

    plot_sol(fig,333,df["x"], df["y"],XY[:,0],XY[:,1],total.real,"$u_{tot}$ : real",cmap)
    plot_sol(fig,336,df["x"], df["y"],XY[:,0],XY[:,1],total.imag,"$u_{tot}$ : imag",cmap)
    plot_sol(fig,339,df["x"], df["y"],XY[:,0],XY[:,1],np.absolute(total),"$u_{tot}$ : abs",cmap)

    plt.tight_layout()
    sol_png = file_dir+"solution-"+suffix+".png"
    plt.savefig(sol_png)
    
    j += 1

    plt.show()


# ## Norm

# In[15]:


#function
x0 = 0.0
y0 = 0.0
a = 1.0
b = 1.0
radius = 1.0
def func(t):
    #print(cpz)
    cpx = a*np.cos(2*np.pi*t)+x0#np.random.rand(N)
    cpy = b*np.sin(2*np.pi*t)+y0#np.random.rand(N)
    out = np.zeros(shape=(len(t),2))
    for i in range(len(t)):
        out[i,0] = cpx[i]
        out[i,1] = cpy[i]
    return out


# In[17]:


#LN
LN = pd.DataFrame(data=Lebesgue_norm,index=N_arr_2,columns=n_arr)
LN.to_csv("files/BEM/circle-convergence/lebesgue_norm.csv",index_label="index")
LN.head()


# In[18]:


#LNv
LNv = pd.DataFrame(data=Lebesgue_norm_var,index=N_arr_2,columns=n_arr)
LNv.to_csv("files/BEM/circle-convergence/lebesgue_norm_var.csv",index_label="index")
LNv.head()


# In[16]:


#calcolo la norma
def get_float(txt):
    return re.findall("[^a-zA-Z:]([-+]?\d+[\.]?\d*)", txt)

#
#i = 0
#n = 4
N_arr_2 = list(np.arange(10,21))+list(np.arange(30,110,10))
n_arr = [4]
Lebesgue_norm = np.zeros((len(N_arr_2),len(n_arr)))
Lebesgue_norm_var = np.zeros((len(N_arr_2),len(n_arr)))

for k in range(len(n_arr)):
    for i in range(0,len(N_arr_2)):    

        N = N_arr_2[i]
        n = n_arr[k]
        count = str(k)+"/"+str(i)+"       "

        print(count,end="\r")

        #
        file_dir = "files/BEM/circle-convergence/n="+str(n)+"/"
        suffix = "P="+str(P)+"-N="+str(N)+"-k="+str(wavevector)+"-n="+str(n)

        df = pd.read_csv(file_dir+"grafico-plane_wave-"+suffix+".csv")
        df.set_index("index",drop=True,inplace=True)
        df["inc"]  = [np.complex(j) for j in df["inc"]]
        df["scat"] = [np.complex(j) for j in df["scat"]]
        df["tot"]  = [np.complex(j) for j in df["tot"]]

        #
        kv = Bs.periodic_kv(xminBs,xmaxBs,p=P,n=N)

        #alloco la Bspline
        bs = Bs.Bspline(sh,[kv],properties={"periodic":[True]})

        t = np.linspace(0,1,N+1,endpoint=True)#[0:-2]
        cp = func(t)
        for j in range(len(t)):
            bs.set_cp(j,cp[j])

        #
        T = np.linspace(xminBs,xmaxBs,1000,endpoint=True)
        xy   = bs.evaluate(T)
        uinc = plane_wave(xy)

        #
        XY = np.asarray([ get_float(j) for j in df["xy"] ]).astype(float)

        analytic_pd,analytic = esFFT.analytic_solution_circle(uinc,XY,wmin=-20,wmax=20,radius=radius,                                                              wavevector=wavevector,opts={"return":"both"})
        
        analytic_pd.to_csv(file_dir+"analytic-solution-plane_wave-"+suffix+".csv",index_label=index)

        df["diff"] = df["scat"] - analytic

        #print(i)
        #Lebesgue_norm[i,k] = np.mean(np.absolute(df["diff"])**2)
        #Lebesgue_norm_var[i,k] = np.var(np.absolute(df["diff"])**2)
        
        LN.at[i,k] = np.mean(np.absolute(df["diff"])**2)
        LNv.at[i,k] = np.var(np.absolute(df["diff"])**2)
        
print("Finished")


# In[ ]:


#LN = pd.DataFrame(data=Lebesgue_norm,index=N_arr_2,columns=n_arr)
LN.to_csv("files/BEM/circle-convergence/lebesgue_norm.csv",index_label="index")
LN.head()


# In[ ]:


#LNv = pd.DataFrame(data=Lebesgue_norm_var,index=N_arr_2,columns=n_arr)
LNv.to_csv("files/BEM/circle-convergence/lebesgue_norm_var.csv",index_label="index")
LNv.head()


# ## Grafico

# In[19]:


#reading
LN = pd.read_csv("files/BEM/circle-convergence/lebesgue_norm.csv")
LN.set_index("index",drop=True,inplace=True)

LNv = pd.read_csv("files/BEM/circle-convergence/lebesgue_norm_var.csv")
LNv.set_index("index",drop=True,inplace=True)


# In[20]:


#grafico
fig = plt.figure(figsize=(15,5))

#
ax = fig.add_subplot(111)
ax.errorbar(LN.index,LN["4"],yerr=LNv["4"],color="blue",label="n=4",linestyle="--",marker=".")
ax.errorbar(LN.index,LN["6"],yerr=LNv["6"],color="green",label="n=6",linestyle="--",marker=".")
ax.errorbar(LN.index,LN["8"],yerr=LNv["8"],color="red",label="n=8",linestyle="--",marker=".")
#ax.scatter(N_arr_2,Lebesgue_norm,color="red",label="calculations")
#plt.xlim(8,41)
plt.grid(True)
plt.xlabel("d.o.f.")
plt.ylabel("Lebesgue norm $L^2$")
plt.legend()
plt.title("Convergence")

plt.show()


# In[ ]:




