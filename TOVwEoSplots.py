from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from EoS import EoSclass
from scipy.interpolate import interp1d
from TOVwEoS import TOVEoS


# Plot of pressure vs rho 
f1 =plt.figure(1)
E = EoSclass
plt.loglog(E.EoSrho,E.EoSP,color='red')
plt.xlabel('log $\epsilon$')
plt.ylabel('log $P$')
plt.title('Pressure vs Density')
plt.legend(['Pressure vs Densty (EoS)'])
# plt.show()

# pressure vs radius and mass vs radius plots

P0 = 0.4
sol = TOVEoS(P0,1e-1)
P = sol[3]
r_arr = sol[2]
m_arr = sol[4]

f2 = plt.figure(2)
plt.plot(r_arr,P,color= 'green')
plt.yscale('log')
plt.ylabel('$P(r)$')
plt.xlabel('$r$')
plt.legend(['$P(r)$ with EoS'])

f3 = plt.figure(3)
plt.plot(r_arr,m_arr,color='blue')
plt.ylabel('$m(r)$')
plt.xlabel('$r$')
plt.legend(['$m(r)$ with EoS'])

plt.show()

