from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from EoS import EoSclass,EoSIntClass
from scipy.interpolate import interp1d
from TOVwEoS import TOVEoS


# Plot of pressure vs rho 
f1 =plt.figure(1)
y = [0,0.01,1,10,10**2,10**3]
E = EoSclass
EI = EoSIntClass
lines = [ '-', '--', '-.', ':','-.','--']
plt.loglog(E.EoSrho,E.EoSP,color='red')

for i in range(len(y)):
    a = EI(y[i])
    plt.loglog(a.rho,a.P,linestyle=lines[i])
plt.legend(['without interactions','y=0','y=0.01','y=1','y=10','y=$10^2$','y=$10^3$'])
plt.xlabel('log $\epsilon$')
plt.ylabel('log $P$')
plt.title('Pressure vs Density')
plt.show()

# pressure vs radius and mass vs radius plots

# P0 = 0.4
# sol = TOVEoS(P0,1e-1)
# P = sol[3]
# r_arr = sol[2]
# m_arr = sol[4]

# f2 = plt.figure(2)
# plt.plot(r_arr,P,color= 'green')
# plt.yscale('log')
# plt.ylabel('$P(r)$')
# plt.xlabel('$r$')
# plt.legend(['$P(r)$ with EoS'])

# f3 = plt.figure(3)
# plt.plot(r_arr,m_arr,color='blue')
# plt.ylabel('$m(r)$')
# plt.xlabel('$r$')
# plt.legend(['$m(r)$ with EoS'])

# plt.show()

