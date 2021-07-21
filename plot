import numpy as np
import matplotlib.pyplot as plt
from AllTOV import AllTOV

eps = np.finfo(float).eps
P0 = 1e5
rho1 = 1e3
tau = 1e-7
sol = AllTOV(P0,rho1, rho2=100, r_c=0.005, order=["incom","incom"])


plt.plot(sol[3], sol[5], color="red")
plt.plot(sol[4], sol[6],color="blue")
print(sol[4][0])
plt.yscale('log')
plt.show()
