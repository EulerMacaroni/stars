import numpy as np
import matplotlib.pyplot as plt
from TOV import TOV
from functions import rel, mass
from TOV2step import TOV2step

eps = np.finfo(float).eps
P = np.linspace(1e-5,1e5,100)
rho = 5000           
r_m = np.array([])
for i in range(len(P)):
    print(i)
    sol = TOV(P[i],rho)
    r_m = np.append(r_m,sol[5])
    if sol[5] <= 9/4:
        print('Illegal')

plt.ylabel('Central Pressure, $P_0$')
plt.xlabel('R/M')
plt.xlim([0,20])
plt.title('Central Pressure vs R/M with rho = %f'%rho)
plt.axvline(x=9/4,label='Buchdahl limit',color='green',linestyle = 'dashed')
plt.plot(r_m,P)
plt.show()



