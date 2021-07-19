from numpy.core.function_base import linspace
from functions import findXPoint
import numpy as np
import matplotlib.pyplot as plt
from TOVwEoS import TOVEoS,rhof
from EoS import EoSclass

E = EoSclass

def Pforrho(rho_0):
    i = min(E.EoSrho,key=lambda x:abs(x-rho_0)) #find closest value to P0
    a = np.where(E.EoSrho==i) # index of closest point to P0
    index =a[0][0] #index of closest P0 (a outputs 2 dim. array)
    
    x2 = E.EoSP[index+1]
    x1 = E.EoSP[index]
    y2 = E.EoSrho[index+1]
    y1 = E.EoSrho[index]

    x3 = findXPoint(x1,x2,y1,y2,rho_0)
    return x3

rho = linspace(10**(-8),10,5,endpoint=True)
P0 = np.array([])

for i in range(len(rho)):
    P0 = np.append(P0,Pforrho(rho[i]))


# P0 = np.linspace(0.9,4,4,endpoint=True)        # limits = (1e-78,4,10)

R = np.array([])
M = np.array([])
for i in range(len(P0)):
    print(i,'with P=',P0[i])
    sol = TOVEoS(P0[i])
    # plt.plot(sol[2],sol[4],color='red')
    plt.plot(sol[0],sol[1],'.',color='black')
    M = np.append(M,sol[1])
    R = np.append(R,sol[0])

plt.xlabel('dimensionaless $r$')
plt.ylabel('dimensionaless $m(r)$')

plt.show()