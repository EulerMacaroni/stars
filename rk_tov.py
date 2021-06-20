# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 15:05:09 2021

@author: Benutzer01
"""

import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import pandas as pd
from scipy.integrate import odeint

eps = np.finfo(float).eps
G= 1
c=1
pi=np.pi
rho_0 = 0.001 #convED(10**5) #density of fluid in geom. units   (10**37)
r = np.arange(eps,14,0.0001)  # time steps   #radius of neutron star 10**3
rmin=0.01
rmax=12. #100

print('M/R limit =', 4/9)
print('rho_0=',rho_0)
print('eps=',eps)

def rho(r):
    rho = rho_0
    return 

def dm(r,Mn):
    dm = 4 * np.pi * rho_0 * r**2 
    return dm

def dP(r,Pn,Mn,G,rho):   #TOV
    dP = -G*Mn*rho_0/r**2 *(1.+Pn/rho_0)*(1.+3.*Pn/rho_0)/(1.-2.*G*Mn/r)
    return dP

    
n=100       
hstep = (rmax-rmin)/float(n)  
r0=rmin
rlist = [r0]
Rstar=0.

mst0=eps
#Mn=0.
masslist = [mst0]
#P0=2.8659734E-03 #0.0028659
s=10
Pst0=5E-4
#Pn=0.
Patmo=1E-8
pressurelist = [Pst0]
massfin = {}   # trying to make a dir and add the lists for different j 

for j in range(0,s):    # loop over P for different P0 to find region of R 
   Pn=0.
   Mn=0.
   m0=mst0
   P0=Pst0*np.real(j)
   #massfin[j]=masslist
   for i in range(0,n):    # loop over r --> calculate runge kutta m and p for every r
     r0=rmin+np.real(i)*hstep
     rlist.append(r0)
     k1 = hstep * (dm(r0, m0))
     k2 = hstep * (dm((r0+hstep/2), (m0+k1/2)))
     k3 = hstep * (dm((r0+hstep/2), (m0+k2/2)))
     k4 = hstep * (dm((r0+hstep), (m0+k3)))
     k = (k1+2*k2+2*k3+k4)/6
     Mn = m0 + k
     m0=Mn
     print('%.4f\t%.4f\t%.4f'% (r0,m0,Mn) )
     print('-------------------------')
     masslist.append(m0)
     k1 = hstep * (dP(r0, P0,Mn,G,rho))
     k2 = hstep * (dP((r0+hstep/2), (P0+k1/2),Mn,G,rho))
     k3 = hstep * (dP((r0+hstep/2), (P0+k2/2),Mn,G,rho))
     k4 = hstep * (dP((r0+hstep), (P0+k3),Mn,G,rho))
     k = (k1+2*k2+2*k3+k4)/6
     Pn = P0 + k
     P0=Pn
     print('%.4f\t%.4f\t%.4f'% (r0,P0,Pn) )
     print('-------------------------')
     pressurelist.append(P0)
     if Pn <= Patmo:
         Rstar=r0
         print('found Rstar',Rstar)
         break

#print('this is the masslist for j=2',massfin[2])
#print(rlist)
#print(masslist)

plt.plot(rlist,masslist)  #massfin[1])     # plotting mass
#plt.plot(rlist,pressurelist)

lss = ['-', '--', '-.', ':']
clr = ['tab:red','tab:blue', 'tab:orange', 'tab:green']
     
fig, ax = plt.subplots()      # plotting pressure

ax.plot(rlist,pressurelist,  color=clr[0])
ax.set_xlabel(r'Radius r [km]')
ax.set_xlim(0,12)

ax.set_ylabel(r'Pressure P')
ax.set_ylim(0,0.004)

ax.set_title(r"The pressure of an incompressible fluid")
ax.legend(title="$Pressure$", frameon=False)

plt.minorticks_on()