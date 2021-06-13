import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions1 import mass,Mass,Nrel
##
rho_0 = 10**1 #density of fluid
G = 6.67408*10**(-11)
r = np.arange(np.finfo(float).eps,1000,10)  # time steps
R= r[len(r)-1] # radius of star 

y0 = 2/3 * np.pi * G * rho_0**2 * R**2  #initial pressure condition at r=0
# y0 = 
print(y0)
# function that returns dy/dt
def fun(P,r):
    k = G*mass(r,rho_0)*rho_0
    dydt = -1*k*(1/(r**2))
    return dydt

# solve ODE
P = odeint(fun,y0,r)

# Uses P equation to find exact values
real_P = []
for item in r:  
    e = y0 * (1 - item**2 / R**2)
    # e = (1 - item**2 / R**2) 
    real_P.append(e)

# calculate error
error = []
for i in range(len(r)):
    error_i = real_P[i]/P[i] -1
    error.append(error_i)

# plot results
plt.plot(r,P)
plt.plot(r,real_P, ':')
plt.plot(r,error)
plt.xlabel('radius $r$')
plt.ylabel('Pressure $P$')
plt.legend(["Calculated Pressure", "Exact Pressure","Error"], loc="upper right")
plt.ticklabel_format(useOffset=False, style='plain')
plt.show()