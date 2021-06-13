import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
##
rho_0 = 10**12#density of fluid
G = 6.67408*10**(-11)
# r = np.linspace(np.finfo(float).eps,10,100)  # time steps
r = np.arange(np.finfo(float).eps,1000,0.1)  # time steps
R= r[len(r)-1] # radius of star 

y0 = 2/3 * np.pi * G * rho_0**2 * R**2  #initial pressure condition at r=0
print(y0)

def mass(r):
    return (4*np.pi/3)*rho_0*(r**3)

# function that returns dy/dt
def model(P,r):
    k = G*mass(r)*rho_0
    dydt = -1*k*(1/(r**2))
    return dydt

# solve ODE
P = odeint(model,y0,r)

#Uses P equation to find exact values
real_P = []
for item in r:  
    e = y0 * (1 - item**2 / R**2)
    real_P.append(e)

# plot results
plt.plot(r,P)
plt.plot(r,real_P, ':')
plt.xlabel('radius')
plt.ylabel('Pressure')
plt.legend(["Calculated Pressure", "Exact Pressure"], loc="upper right")
plt.ticklabel_format(useOffset=False, style='plain')
plt.show()