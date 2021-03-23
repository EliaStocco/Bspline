#!/usr/bin/env python
# coding: utf-8

# In[8]:


get_ipython().system('jupyter nbconvert --to script FFT.ipynb')


# # FFT

# In[50]:


from scipy.fft import fft, ifft, fftfreq, fftshift
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def FFT(y,x=None,opts=None):
    
    if opts is None:
        opts = {}
    if "xmin" not in opts:
        opts["xmin"] = 0.
    if "xmax" not in opts:
        opts["xmax"] = 1.
    if "plot" not in opts:
        opts["plot"] = False
    if "fig" not in opts:
        opts["fig"] = None
    if "ax0" not in opts or "ax1" not in opts:
        opts["ax0"] = 121
        opts["ax1"] = 122
    if "inv" not in opts:
        opts["inv"] = False
    
    y = np.asarray(y)
    N = len(y)
    
    if opts["inv"] == False :
        yf = fft(y,norm="ortho")
        xf = fftshift(fftfreq(N)*N).astype(int)
        yplot = fftshift(yf)   
    else :
        yf = fftshift(y)   
        yplot = ifft(yf,norm="ortho")
        xf = fftshift(fftfreq(N)*N).astype(int)
        
     
    out = pd.DataFrame(data=yplot,index=xf,columns=["fft"],dtype=np.complex)
    
    if opts["plot"] == True:
        
        if x is None:
            x = np.linspace(xmin,xmax , N, endpoint=False)
        #
        if opts["fig"] is None:
            fig = plt.figure ( 0 , figsize = ( 15 , 5 ) )
        else :
            fig = opts["fig"]
            
        if opts["inv"] == False:
            x1 = x
            x2 = xf
            title1 = "Original Function"
            title2 = "Fourier Transform"
        else :
            x1 = xf
            x2 = x
            title1 = "Fourier Transform"
            title2 = "Original Function (Anti-Trasform)"

        #
        ax = fig.add_subplot(opts["ax0"])
        plt.plot(x1, y.real,color="blue" ,label="real")#,marker="+")
        plt.plot(x1, y.imag,color="green",label="imag")#,marker="x")
        #plt.xlim(min(x),max(x))
        plt.legend()
        plt.grid(True)
        plt.title(title1)

        #
        ax = fig.add_subplot(opts["ax1"])
        plt.plot(x2, yplot.real,color="blue" ,label="real")#,marker="+")
        plt.plot(x2, yplot.imag,color="green",label="imag")#,marker="x")
        #plt.xlim(min(xf),max(xf))
        plt.legend()
        plt.grid(True)
        plt.title(title2)
        
        plt.show()

    return out


# In[51]:


import numpy as np

xmin = 0.
xmax = 1.
N=1000

x = np.linspace(xmin,xmax , N, endpoint=False)
w = 10
y = np.exp(1.j * w * x *2*np.pi)# + np.exp(1.j * 3 * x * 2*np.pi)

out = FFT(y,x,opts={"plot":True})


# In[54]:


#inv = np.asarray(out["fft"])
#xinv = np.asarray(out.index)
orig = FFT(out["fft"],x,opts={"plot":True,"inv":True})


# In[128]:


import scipy.special

wavevector = 30
radius = 1.0
k = wavevector
R = radius
def H_coeff(l,r):
    return scipy.special.hankel1(l,k*r)/scipy.special.hankel1(l,k*R)

r_array = np.linspace(4,10,100)

analytic_x = pd.DataFrame(index=r_array,columns=out.index,dtype=np.complex)

for r in r_array:
    for i in out.index:
        #print(i,"-",r,end="\r")
        analytic_x.at[r,i] = - out.at[i,"fft"] * H_coeff(i,r)


# In[148]:


l = i
scipy.special.hankel1(499,100)


# In[144]:


k*R


# In[129]:


analytic_x


# In[126]:



analytic = pd.DataFrame(index=r_array,columns=out.index,dtype=np.complex)
analytic_np = np.asarray(analytic)
analytic_x_np = np.asarray(analytic_x)
for r in range(len(r_array)):
    analytic_np[r,:] = FFT(analytic_x_np[r,:],opts={"inv":True})["fft"]
analytic = pd.DataFrame(data=analytic_np,index=r_array,columns=out.index,dtype=np.complex)


# In[127]:


analytic


# In[ ]:




