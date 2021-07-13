import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions import rel, mass, convED,convP
from TOV import TOV

eps = np.finfo(float).eps

rho = 1e3
P0 = 10e8
tau = 1e-5
sol = TOV(P0,rho,tau)

R = sol[0]
M = sol[1]
r_arr = sol[2]
m_arr = sol[4]
P = sol[3]
M2 = mass(R,rho)


ExactP =  rel(r_arr,R,M,rho)
ExactM = mass(r_arr,rho)
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Radius $r$')
ax1.set_ylabel('Pressure $P$', color=color)
plt.yscale('log')
ax1.plot(r_arr, P, color='green')
ax1.plot(r_arr, ExactP, color=color,linestyle = 'dashed')
ax1.legend(['Numerical P(r)','Analytical P(r)'])

ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Mass $M$', color=color)  # we already handled the x-label with ax1
ax2.plot(r_arr, m_arr, color=color)
ax2.plot(r_arr, ExactM, color='darkviolet',linestyle = 'dashed')
ax2.legend(['Numerical m(r)','Analytical m(r)'])


ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()