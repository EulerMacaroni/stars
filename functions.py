import numpy as np
pi = np.pi

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


# funtion for equation of state 

def EosRho(k):
    m_f =1
    return (1/(pi**2))*(k**2)*(np.sqrt(m_f**2 + k**2))
    
def EosP(k):
    m_f =1
    return (1/(3*pi**2))*((k**4)/(np.sqrt(m_f**2 + k**2)))

def intTerm(y,z):
    return (((1/(3*pi**2)))**2)*(y**2)*(z**6)

def findXPoint(x1,x2,y1,y2,y3):
    m = (x2 - x1)/(y2 - y1)
    x3 = (y3 - y1)*m + x1
    return x3


def tauf1(P0):
    if P0 <= 1e-7:
        return 2
    elif P0> 1e-7 and P0 <= 1e-6:
        return 0.5
    elif P0> 1e-6 and P0 <= 1e-5:
        return 0.25
    elif P0> 1e-5 and P0 <= 1e-4:
        return 0.1
    elif P0 > 1e-4 and P0 <= 1e-3:
        return 0.01
    elif P0 > 1e-3 and P0 <= 1e-2:
        return 0.005
    elif P0 > 1e-2 and P0 <= 1e-0:
        return 0.0005
    elif P0 > 1e-0 and P0 <= 1e2:
        return 0.0001
    elif P0 > 1e2 and P0 <= 1e3:
        return 0.00001
    else:
        return P0*(1e-8)



# def tauf2(rho):
#     if rho < 1:
#         return 0.1/2
#     elif rho >= 1:
#         return 0.1

def tauf2(rho):
    return 0.00001
