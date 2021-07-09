import numpy as np
from numpy.core.function_base import linspace
import matplotlib.pyplot as plt
from EoSTOV import TOVEOS

m_r_array = linspace(1e-8, 4/9 - 1e-5, num = 25)
rho = 1
R_array = []
M_array = []
for m_r in m_r_array:
    sol = TOVEOS(rho, m_r)
    R_array.append(sol[0])
    M_array.append(sol[1])

R = sol[0]
M = sol[1]
#r_arr = sol[2]



fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Radius $R$')
ax1.set_ylabel('Mass $M$', color=color)
# plt.yscale('log')
ax1.plot(R_array, M_array, color='green')

ax1.tick_params(axis='y', labelcolor=color)

plt.show()
