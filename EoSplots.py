from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from EoS import EoSclass,EoSIntClass
from scipy.interpolate import interp1d
from TOVwEoS import TOVEoS


# Plot of pressure vs rho 

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