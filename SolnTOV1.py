# from TOV import M_R
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from functions1 import convED ,mass , rel , Nrel, Mass

eps = np.finfo(float).eps
G= 1
c=1
rho_0 = convED(10**5) #density of fluid in geom. units   (10**37)
r = np.arange(eps,100,0.01)  # time steps   #radius of neutron star 10**3
m = mass(r,rho_0)
m_r = np.divide(m,r)
R = r[len(r)-1]
M = Mass(R,rho_0)
M_R = M/R          # M/R < 4/9
print('mass=',M)
print('radius=', R)
print('M/R=',M_R,'M/R limit =', 4/9)
print('rho_0=',rho_0)
print('eps=',eps)

def fun(r, P):
    k = G*mass(r,rho_0)*rho_0
    factor_1 = 1 + P/(rho_0 * c**2)
    factor_2 = (1 + (4 * np.pi * r**3 * P)/(mass(r,rho_0) * c**2))
    factor_3 = 1 / (1 - (2 * G * mass(r,rho_0))/(r * c**2))
    dPdr = -1*k*(1/(r**2)) * factor_1 * factor_2 * factor_3
    return dPdr

tspan = np.arange(eps,100,0.01)
# yinit = [0]
yinit = [rho_0*(c**2)*(np.sqrt(1-M_R)-1)/(1-3*np.sqrt(1-M_R))]
# c1= []
# sol = solve_ivp(fun, tspan , yinit, args=(1.5, 1, 3, 1),
                # dense_output=True)

sol = solve_ivp(fun, r,yinit, 
                dense_output=True)
t = r
z = sol.sol(t)
plt.plot(t, z.T)
# plt.plot(sol.t,sol.y[0])
plt.xlabel('radius')
plt.ylabel('pressure')
# plt.plot(t,Nrel(yinit,r,R))
# plt.plot(t,rel(r,R,M,rho_0))
plt.legend(['Num', 'Rel'], shadow=True)
# plt.title('Lotka-Volterra System')
plt.show()