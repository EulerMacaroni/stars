import numpy as np
import matplotlib.pyplot as plt
from TOV import TOV
from functions import rel, mass
from TOV2step import TOVstep

eps = np.finfo(float).eps
P = np.linspace(1e1,1e15,10)
rho = 5000         
R1   = np.array([])
R2   = np.array([])
M1   = np.array([])
M2   = np.array([])
r_m1 = np.array([])
r_m2 = np.array([])

for i in range(len(P)):
    print(i)
    sol2 = TOVstep(P[i],rho)
    sol1 = TOV(P[i],rho)
    r_m1 = np.append(r_m1,sol1[5])
    r_m2 = np.append(r_m2,sol2[5])
    R1   = np.append(R1,sol1[0])
    R2   = np.append(R2,sol2[0])
    M1   = np.append(M1,sol1[1])
    M2   = np.append(M2,sol2[1])

    if sol1[5] <= 9/4:
        print('Illegal for 1 step profile')

plt.ylabel('Central Pressure, $P_0$')
plt.xlabel('R/M')
plt.xlim([0,5])
plt.title('Central Pressure vs R/M with rho = %f'%rho)
plt.axvline(x=9/4,label='Buchdahl limit',color='green',linestyle = 'dashed')
plt.plot(r_m1,P,color='red')
plt.plot(r_m2,P,color='black')
plt.legend(['Buchdahl limit','1 Step Profile','2 Step Profile'])
plt.show()

diff = np.subtract(r_m2,r_m1)
print('differnce b/w r_m2 and r_m1',diff)
# plt.plot(diff,'.',color='darkviolet')
# plt.show()

# plt.plot(R1,M1)
# plt.plot(R2,M2)
# plt.xlabel('R')
# plt.ylabel('M')
# plt.legend(['R vs M (1 Step)','R vs M (2 Step)'])
# plt.show()

plt.plot(R1,np.power(r_m1,-1),'.',color= 'green')
plt.plot(R2,np.power(r_m2,-1),'.',color='red')
plt.axhline(y=4/9,color='green',linestyle = 'dashed')
plt.legend(['1 Step','2 Step','M/R < 4/9'])
plt.xlabel('R')
plt.ylabel('M/R')
plt.ylim([0.2,4/9 +0.1 ])
plt.grid()
plt.show()
