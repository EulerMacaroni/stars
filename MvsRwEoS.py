from functions import EosP
import numpy as np
import matplotlib.pyplot as plt
from TOVwEoS import TOVEoS
from EoS import EoSclass

P0 = np.linspace(1e-7,1,10,endpoint=True)        # limits = (1e-78,4,10)
def tauf(P):
    if P < 0.5:
        return 1e-1
    elif P>=0.5:
        return 1e-3

R = np.array([])
M = np.array([])
for i in range(len(P0)):
    print(P0[i])
    sol = TOVEoS(P0[i],tauf(P0[i]))
    M = np.append(M,sol[1])
    R = np.append(R,sol[0])

plt.plot(R,M,color='black')
plt.xlabel('dimensionaless R')
plt.ylabel('dimensionaless M')
plt.legend(['$R$ vs $M$'])
plt.show()
# E = EoSclass

# rho= 10
# i = min(E.EoSrho,key=lambda x:abs(x-rho)) #find closest value to P0
# a = np.where(E.EoSrho==i) # index of closest point to P0
# print(a)
# index =a[0][0] #index of closest P0 (a outputs 2 dim. array)
# print('for given rho = ',rho,'Pressure =', E.EoSP[index])