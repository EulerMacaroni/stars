import numpy as np
import matplotlib.pyplot as plt
from EoS import EoSclass
from scipy.interpolate import interp1d
from TOVwEoS import TOVEoS

# pressure vs radius and mass vs radius plots

P0 = 1e-4
sol = TOVEoS(P0,1e-2)
P = sol[3]
r_arr = sol[2]
m_arr = sol[4]

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Radius $r$')
ax1.set_ylabel('Pressure $P$', color=color)
plt.yscale('log')
ax1.plot(r_arr, P, color='green')
ax1.legend(['Numerical P(r)'])

ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Mass $M$', color=color)  # we already handled the x-label with ax1
ax2.plot(r_arr, m_arr, color=color)
ax2.legend(['Numerical m(r)'])


ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
