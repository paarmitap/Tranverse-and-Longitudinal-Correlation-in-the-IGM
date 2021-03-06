#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 09:42:01 2019

@author: esha
"""


import numpy as np 
import matplotlib.pyplot as plt
import glob


#list_B=sorted(glob.glob("Pair_data/*/*B.dat"))

#list_A=sorted(glob.glob("Pair_data/*/*A.dat"))

#theta=np.genfromtxt("all_pairs_address_file.txt",unpack=True,usecols=11)
#tht1=np.genfromtxt("all_pairs_address_file.txt",unpack=True,usecols=11)

#list_r=sorted(glob.glob("Pair_data/*/*bg*"))
#list_r1=sorted(glob.glob("Pair_data/*/*fg*"))


list_B=sorted(glob.glob("*B*"))

list_A=sorted(glob.glob("*A*"))

theta=np.genfromtxt("all_pairs_address_file.txt",unpack=True,usecols=11)
tht1=np.genfromtxt("all_pairs_address_file.txt",unpack=True,usecols=11)

list_r=sorted(glob.glob("*bg*"))
list_r1=sorted(glob.glob("*fg*"))

z1,z2=np.genfromtxt("all_pairs_address_file.txt",unpack=True,usecols=(4,9))





total=len(list_A)

from scipy.interpolate import CubicSpline
g=np.zeros(total,)


for i in range(total):
    q1=[]
    l1=[]
    n1=[]
    
    q1=np.loadtxt(list_A[i])
    
    
    plt.plot(q1[:,0],q1[:,1])
    l1=q1[:,0]
    n1, = np.where((l1<=(1215.6701-12.16)*(1+z2[i])) & (l1>=(1025.7223+3.416)*(1+z1[i])))
    reg1=[]
    reg1=np.take(q1,n1,axis=0)
    plt.plot(reg1[:,0],reg1[:,1],color='r')
    
    a1,a2= np.genfromtxt(list_r[i],unpack=True,usecols=(0,1))

    b1=[]
    b1=reg1[:,0]

    dele1=[]
    for j in range(len(a1)):   

        for k in range(len(b1)):
            if((b1[k]<=a2[j]) & (b1[k]>=a1[j])):
                    dele1.append(k)
    NEWreg1=np.delete(reg1,dele1,axis=0) 
               
    lam1=[]
    lam1=NEWreg1[:,0]
    
    zA1=[]
    zA1=(lam1/1215.67)-1
    
    # converting zabs to velocity
    vel1=[]
    vel1= ((((1+z2[i])**2)-(1+zA1)**2)/(((1+z2[i])**2)+(1+zA1)**2))*299792
    
    # finding normalised flux
    normal1=[]
    normal1=NEWreg1[:,1]/NEWreg1[:,2]
    m1=np.mean(normal1)
    normal1=normal1-m1
    
    

    q2=[]
    q2=np.loadtxt(list_B[i])
    l2=[]
    l2=q2[:,0]
    n2=[]
    n2, = np.where((l2<=(1215.6701-12.16)*(1+z2[i])) & (l2>=(1025.7223+3.416)*(1+z1[i])))
    reg2=[]
    reg2=np.take(q2,n2,axis=0)
    plt.plot(q2[:,0],q2[:,1],color='y')
    plt.plot(reg2[:,0],reg2[:,1],color='g')
    plt.legend([list_B[i]],fontsize=12)
    plt.show()
    a3,a4= np.genfromtxt(list_r1[i],unpack=True,usecols=(0,1))
    b2=[]
    b2=reg2[:,0]

    dele2=[]
    for j in range(len(a3)):   

        for k in range(len(b2)):
            if((b2[k]<=a4[j]) & (b2[k]>=a3[j])):
                    dele2.append(k)
    NEWreg2=np.delete(reg2,dele2,axis=0) 
               
    lam2=[]
    lam2=NEWreg2[:,0]
    
    zA2=[]
    zA2=(lam2/1215.67)-1
    
    # converting zabs to velocity
    vel2=[]
    vel2= ((((1+z2[i])**2)-(1+zA2)**2)/(((1+z2[i])**2)+(1+zA2)**2))*299792
    
    # finding normalised flux
    normal2=[]
    normal2=NEWreg2[:,1]/NEWreg2[:,2]
    m2=np.mean(normal2)
    normal2=normal2-m2
    
    
    vel2= vel2[::-1] #np.flip(vel2,axis=None)
    normal2=normal2[::-1] #np.flip(normal2,axis=None)
    vel1=vel1[::-1] #np.flip(vel1,axis=None)
    normal1=normal1[::-1] #np.flip(normal1,axis=None)
   
    cs1=[]
    cs1 = CubicSpline(vel2, normal2)
    b=cs1(vel1)
    g[i]=np.corrcoef(normal1,b)[0,1]
    
    plt.plot(vel1,normal1,'o',markersize=1,color='r')
    plt.plot(vel2,normal2,'o',markersize=2,color='b')
    plt.plot(vel1,b,'o',markersize=1,color='g')
    plt.xlabel("vel")
    plt.ylabel("Flux")
    plt.legend([list_B[i]],fontsize=12)
    plt.show()


# Finding angular diameter distance  
r=np.zeros(total,)
vperp=np.zeros(total,)
import scipy.integrate
f=lambda x: 1/(np.sqrt(0.3*(1+x)**3 + 0.7))

for i in range(total):
    r[i],_=scipy.integrate.quad(f, 0, z2[i])
    
for i in range(total):
    r[i]=(r[i]*3000)/(1+z2[i])
    

theta=theta*0.00029

for i in range(total):
    vperp[i]=r[i]*100*(np.sqrt(0.3*(1+z2[i])**3 + 0.7))*theta[i]



#result=np.zeros((total,3))
#for i in range(total):
#    result[i,0]=tht1[i]
#    result[i,1]=g[i]
#    result[i,2]=vperp[i]
#np.savetxt('result2.txt',result,delimiter='\t')

'''
import scipy
from scipy import stats


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

a1,b1,c1,d1=split(g,4)
a2,b2,c2,d2=split(vperp,4)
gbin=np.array([np.median(a1),np.median(b1),np.median(c1),np.median(d1)])
vbin=np.array([np.median(a2),np.median(b2),np.median(c2),np.median(d2)])
plt.plot(vbin,gbin,'o')

'''



fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()



ax1.plot(vperp, g,'o')
ax2.plot(tht1,g,'o') # Create a dummy plot
ax2.lines=[]
ax1.set_xlabel('r (km/s)')
ax1.set_ylabel('corr coeff')
ax2.set_xlabel('Theta in arcmin')

plt.show()