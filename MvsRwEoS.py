from functions import EosP
import numpy as np
import matplotlib.pyplot as plt
from TOVwEoS import TOVEoS
from EoS import EoSclass

P0 = np.linspace(1e-7,10**1,50,endpoint=True)        # limits = (1e-78,4,10)
# tau =1e-3
R = np.array([])
M = np.array([])
for i in range(len(P0)):
    print(P0[i])
    sol = TOVEoS(P0[i])
    plt.plot(sol[2],sol[4],color='red')
    plt.plot(sol[0],sol[1],'.',color='black')
    M = np.append(M,sol[1])
    R = np.append(R,sol[0])

plt.xlabel('dimensionaless $r$')
plt.ylabel('dimensionaless $m(r)$')

plt.show()