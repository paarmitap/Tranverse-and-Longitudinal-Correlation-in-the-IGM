#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:43:40 2019

@author: esha
"""

#from pylab import *
import numpy as np
import matplotlib.pyplot as plt



#Getting values of z
z=np.genfromtxt("all_pairs_address_file.txt",usecols=9)

# Getting list of all B quasars
import glob
list_B=sorted(glob.glob("*B*"))

for i in range(29):
    q1=[]
    q1=np.loadtxt(list_B[i])
    l1=[]
    l1=q1[:,0]
    n2=[]
    n2, = np.where((l1<=1215.6701*(1+z[i])) & (l1>=1025.7223*(1+z[i])))
    
    v1=[]
    v1=np.take(q1,n2,axis=0)
    np.savetxt(list_B[i],v1)


#Converting wavlength into array of Zabs
for i in range(29):
    r=[]
    r=np.loadtxt(list_B[i])
    s=[]
    s=r[:,0]
    zA=[]
    zA=(s/1215.67)-1
    c=[]
    c=np.delete(r,0,axis=1)
    t=[]
    t=np.insert(c,0,zA,axis=1)
    np.savetxt(list_B[i],t)  
    


# Converting Zabs to velocity
for i in range(29):
    y1=[]
    y1=np.loadtxt(list_B[i])
    w=[]
    w=y1[:,0]
    vel= ((((1+z[i])**2)-(1+w)**2)/(((1+z[i])**2)+(1+w)**2))*299792
    v1=[]
    m1=[]
    v1=np.delete(y1,0,axis=1)
    m1=np.insert(v1,0,vel,axis=1)
    np.savetxt(list_B[i],m1)   
    

# Finding normalized flux
for i in range(29):
    p=[]
    p=np.loadtxt(list_B[i])
    
    normal=[]
    normal=p[:,1]/p[:,2]
    m=np.mean(normal)
    normal=normal-m
    n=[]
    n=np.insert(p,4,normal,axis=1)
    np.savetxt(list_B[i],n)  

#Interpolation and correlation 
from scipy.interpolate import CubicSpline
arr=np.arange(0,20000,100)
lag=range(100)
g=np.zeros((200,29),dtype=float)

for i in range(29):
    a=np.zeros((200,))
    u=[]
    vel=[]
    u=np.loadtxt(list_B[i])
    vel=u[:,0]
    vel=np.flip(vel,axis=None)
    normal=[]
    normal=u[:,4]
    normal=np.flip(normal,axis=None)
    cs = CubicSpline(vel, normal)
    a=cs(arr)

    
    for n in range(100):
             g[n,i]=np.sum((a[0:200-lag[n]])*(a[lag[n]:200]))
             var=np.sum(a**2)
             g[n,i]=g[n,i]/var

        
# Finding median of correlations         

new=np.zeros(100,)
for i in range(100):
    new[i]=np.median(g[i,:],axis=0)
 

plt.plot(range(100),new)
plt.xscale('log')
plt.yscale('log')
plt.show()






