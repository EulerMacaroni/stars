from numpy.core.function_base import linspace
from functions import findXPoint
import numpy as np
import matplotlib.pyplot as plt
from TOVwEoS import TOVEoS,rhof
from EoS import EoSclass
from scipy.interpolate import interp1d

E = EoSclass

def Pforrho(rho):
    i = min(E.EoSrho,key=lambda x:abs(x-rho)) #find closest value to P0
    a = np.where(E.EoSrho==i) # index of closest point to P0
    index =a[0][0] #index of closest P0 (a outputs 2 dim. array)
    f = interp1d([E.EoSrho[index-1],E.EoSrho[index],E.EoSrho[index+1]],[E.EoSP[index-1],E.EoSP[index],E.EoSP[index+1]])
    return f(rho)

rho = np.linspace(10**(-8),10,5,endpoint=True)
P0 = np.array([])

for i in range(len(rho)):
    P0 = np.append(P0,Pforrho(rho[i]))

print(rho)
print(P0)
# P0 = np.linspace(0.9,4,4,endpoint=True)        # limits = (1e-78,4,10)

R = np.array([])
M = np.array([])
for i in range(len(P0)):
    print('Star #',i,'with P=',P0[i],'and rho=',rho[i])
    sol = TOVEoS(P0[i])
    # plt.plot(sol[2],sol[4],color='red')
    M = np.append(M,sol[1])
    R = np.append(R,sol[0])
plt.xlim([0,20])
plt.ylim([-0.05,0.5])

plt.plot(R,M,'.',color='black')
plt.plot(R,M,color='red')

plt.xlabel('dimensionaless $r$')
plt.ylabel('dimensionaless $m(r)$')

plt.show()