import numpy as np
import matplotlib.pyplot as plt
from EoS import EoSclass,EoSIntClass
from scipy.interpolate import interp1d
from TOVwEoS import TOVEoS

# pressure vs radius and mass vs radius plots

P0 = 1
sol = TOVEoS(P0)
P = sol[3]
r_arr = sol[2]
m_arr = sol[4]

f1 = plt.figure(1)
plt.plot(r_arr,P,color= 'green')
plt.yscale('log')
plt.ylabel('$P(r)$')
plt.xlabel('$r$')
plt.legend(['$P(r)$ with EoS'])

f2 = plt.figure(2)
plt.plot(r_arr,m_arr,color='blue')
plt.ylabel('$m(r)$')
plt.xlabel('$r$')
plt.legend(['$m(r)$ with EoS'])

plt.show()
