from interPwEoS import Pforrho
import numpy as np
import matplotlib.pyplot as plt
from TOVwEoS import TOVEoS
from EoS import EoSclass

E = EoSclass

rho = np.linspace(10**(-8),10,5,endpoint=True)
P0 = np.array([])

for i in range(len(rho)):
    P0 = np.append(P0,Pforrho(rho[i]))

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