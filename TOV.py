import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions1 import mass, Mass,Nrel,convED, rel,convM , dP_dr

# define variables
eps = np.finfo(float).eps
rho_0 = convED(10**30) #density of fluid in geom. units   (10**37)
# rho_0 = 10**(-10)
G = 1
c = 1
# solar mass in Geom. units 
M_s = convM(2*(10**40))

# m_r = np.arange(np.finfo(float).eps,1000,0.1)  # time steps
# r = np.linspace(eps,10**3,5000)  # time steps
r = np.arange(eps,100,0.01)  # time steps   #radius of neutron star 10**3
m = mass(r,rho_0)
m_r = np.divide(m,r)
R = r[len(r)-1]
M = Mass(R,rho_0)
# M = ((4/9)*R)-0.1
M_R = M/R          # M/R < 4/9
print('mass=',M)
print('radius=', R)
print('M/R=',M_R,'M/R limit =', 4/9)
print('rho_0=',rho_0)
print('eps=',eps)
# print('array of r/m=', m_r)
# y0 = [eps]
# M_R = 2*Mass(R,rho_0)/R
# M_R = M/R
# print(M_R) 
# print(M)
# M_R  = r/mass(r,rho_0)
# print(M_R) 
y0 = [rho_0*(c**2)*(np.sqrt(1-M_R)-1)/(1-3*np.sqrt(1-M_R))]
# y0 =  [((4/9)-0.1)]
# y0 = [1.8]

print('y0=',y0[0])


# function that returns dP/dt

# def fun(P,r):
#     k = G*mass(r,rho_0)*rho_0
#     factor_1 = 1 + P/(rho_0 * c**2)
#     factor_2 = (1 + (4 * np.pi * r**3 * P)/(mass(r,rho_0) * c**2))
#     factor_3 = 1 / (1 - (2 * G * mass(r,rho_0))/(r * c**2))
#     dPdr = -1*k*(1/(r**2)) * factor_1 * factor_2 * factor_3
#     return dPdr
def dP_dr(P,r):
    k  = 1/(r**2)
    k1 = (rho_0 + (P/(1**2)))
    k2 = (mass(r,rho_0) + (4*np.pi*(r**3)*P)/(1**2))
    k3 = (1 - (2*1*mass(r,rho_0))/(r*(1**2)))**(-1)
    dPdr = -(k)*(k1*k2*k3)
    return dPdr

for i in range(len(y0)):
    P = odeint(dP_dr,y0[i],r)
    plt.plot(r,P)
    # plt.plot(m_r,P)
    # plt.plot(r,real_P)
# plt.plot(r,Nrel(y0,r,R))
# plt.plot(m_r,rel(r,R,M,rho_0))
# plt.legend(["$y_{01}$", "$y_{01}$","$y_{03}$"], loc="upper right")
# plt.legend(["$P_{num}$","$P_{exact}$"])
# plt.legend(["relativistic","non relevistic" ])
# plt.ticklabel_format(useOffset=False, style='plain')
# plt.legend([mass = M])
# plt.axvline(x=4/9)
plt.xlabel('Radius $r$')
plt.ylabel('Pressure $P$')
plt.show()