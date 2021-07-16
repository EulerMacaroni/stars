import numpy as np
import matplotlib.pyplot as plt
from EoS import EoSclass

# Plot of pressure vs rho 
E = EoSclass
plt.loglog(E.EoSrho,E.EoSP,color='red')
plt.xlabel('log $\epsilon$')
plt.ylabel('log $P$')
plt.title('Pressure vs Density')
plt.legend(['Pressure vs Densty (EoS)'])

# P0 = 1e7
# i = min(E.EoSP,key=lambda x:abs(x-P0))
# a = np.where(E.EoSP==i)
# plt.plot(E.EoSrho[a],E.EoSP[a],'o')
# plt.show()