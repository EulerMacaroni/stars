import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions import rel, mass

# single density profile TOV solver

def TOV(P,rho):

    def dm_dr(r,M):
        return (4)*(np.pi)*(r**2)*rho_0

    def dP_dr(r,P):
        k  = (mstep*rho_0)/(r**2)
        k1 = 1 + (P/rho_0)
        k2 = 1 + ((4*np.pi*(r**3)*P)/(mstep))
        k3 = (1 - ((2*mstep)/(r)))**(-1)
        dPdr = -(k)*(k1*k2*k3)
        return dPdr

    eps = np.finfo(float).eps

    rho_0 = rho 
    P0 = P
    
    r_new = 1e-6
    m0 = eps
    t_span =np.linspace(r_new,1e-5+r_new,10)
    
    # # solving Mass and Pressure differential equations

    P_array  = np.array([])
    r_array  = np.array([])
    m_array  = np.array([])
    m_tarray = np.array([])
    ex_array = np.array([])

    while True:
    
        m = odeint(dm_dr,m0,t_span,tfirst=True)
        
        mstep = m[-1]

        Pressure = odeint(dP_dr,P0,t_span,tfirst=True)
        
        P_array = np.append(P_array,Pressure[-1])
        r_array = np.append(r_array,t_span[-1])
        m_array = np.append(m_array,m[-1])
        
        if Pressure[-1] <= 1e-10:
            R = r_array[-1]
            M = m_array[-1]
            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M > 2.25) = ', (R)/(M))
            break
        
        t_span = np.linspace(r_array[-1],1e-5+r_array[-1],10)
    
        P0 = Pressure[-1]
        m0 = m[-1]
    R = r_array[-1]
    M = m[-1]
    compactness = R/M
    return R,M,r_array,P_array,m_array,compactness

