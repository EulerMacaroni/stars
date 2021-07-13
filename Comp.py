import numpy as np
import matplotlib.pyplot as plt
from TOV import TOV
from functions import rel, mass
from TOV2step import TOV2step

eps = np.finfo(float).eps
P0 = np.linspace(1e2,1e7,100)

rho = 1e3       
tau = 1e-5

R1   = np.array([])
R2   = np.array([])
M1   = np.array([])
M2   = np.array([])
r_m1 = np.array([])
r_m2 = np.array([])

for i in range(len(P0)):
    print(i,'with P0=',P0[i])
    sol1 = TOV(P0[i],rho,tau)
    sol2 = TOV2step(P0[i],rho,tau)
    r_m2 = np.append(r_m2,sol2[5])
    r_m1 = np.append(r_m1,sol1[5])
    R2   = np.append(R2,sol2[0])
    R1   = np.append(R1,sol1[0])
    M1   = np.append(M1,sol1[1])
    M2   = np.append(M2,sol2[1])


plt.ylabel('Central Pressure, $P_0$')
plt.xlabel('R/M')
plt.xlim([0,5])
plt.title('Central Pressure vs R/M with rho = %f'%rho)
plt.axvline(x=9/4,label='Buchdahl limit',color='green',linestyle = 'dashed')
plt.plot(r_m1,P0,color='red')
plt.plot(r_m2,P0,color='black')
plt.legend(['Buchdahl limit','1 Step Profile','2 Step Profile'])
plt.show()


plt.plot(R1,np.power(r_m1,-1),'.',color= 'red')
plt.plot(R2,np.power(r_m2,-1),'.',color= 'green')
plt.axhline(y=4/9,color='green',linestyle = 'dashed')
plt.legend(['1 Step','2 Step','Buchdahl limit'])
plt.xlabel('R')
plt.ylabel('M/R')
plt.ylim([0,4/9 +0.1 ])
plt.grid()
plt.show()