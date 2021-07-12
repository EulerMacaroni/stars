import numpy as np

# mass function
def mass(r,rho_0):
    return (4/3)*(np.pi)*(r**3)*(rho_0)

# Mass Function
def Mass(R,rho_0):
    return mass(R,rho_0)

# Non relativistc solution
def Nrel(y0,r,R):
     return y0*(1 - r**2/R**2)

# relativistic solution
def rel(r,R,M,rho_0):
    A = (R*np.sqrt(R-2*M)-np.sqrt((R**3)-2*M*(r**2)))
    B = (np.sqrt((R**3)-2*M*(r**2))-3*R*np.sqrt(R-2*M)) 
    return (rho_0)*(A/B)

def y0(m_r,rho_0):
    
    return rho_0*((np.sqrt(1-2*m_r) -1)/(1-3*np.sqrt(1-2*m_r)))


# converstion energy density (SI --> Geom.)

def convED(rho_0):
    return (rho_0*(6.6743*(10**-11)))/((2.9979*(10**8))**4)

# converstion mass (SI --> Geom.)
def convM(M):
    return (M*(6.6743*(10**-11)))/((2.9979*(10**8))**2)

def convP(P):
    return (P*(6.6743*(10**-11)))/((2.9979*(10**8))**4)