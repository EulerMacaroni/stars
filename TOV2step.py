import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions import rel, mass

# single density profile TOV solver

def TOV2step(P,rho,r_c):

    def rhof(r):
        if r <= r_c:
            return rho
        else: 
            return rho/2 

    def dm_dr(r,M):
        return (4)*(np.pi)*(r**2)*rhof(r)

    def dP_dr(r,P):
        k  = (mstep*rho(r))/(r**2)
        k1 = 1 + (P/rho(r))
        k2 = 1 + ((4*np.pi*(r**3)*P)/(mstep))
        k3 = (1 - ((2*mstep)/(r)))**(-1)
        dPdr = -(k)*(k1*k2*k3)
        return dPdr

    eps = np.finfo(float).eps

    rho_0 = rho 
    P0 = P
    
    r_new = 1e-6 
    m0 = eps
    
    # # solving Mass and Pressure differential equations

    P_array  = np.array([])
    r_array  = np.array([])
    m_array  = np.array([])
    m_tarray = np.array([])
    ex_array = np.array([])

    while True:
    
        t_span = [r_new,(1.001)*r_new]
        # t_span = [r_new+eps]
        
        m = odeint(dm_dr,m0,t_span,tfirst=True)
        mstep = m[-1]
        Pressure = odeint(dP_dr,P0,t_span,tfirst=True)
        P_array = np.append(P_array,Pressure[-1])
        r_array = np.append(r_array,t_span[-1])
        m_array = np.append(m_array,m[-1])
        # print(Pressure[-1])

        if Pressure[-1] <= 1e-10:
            R = r_array[-1]
            M = m_array[-1]
            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M > 2.25) = ', (R)/(M))
            break

        r_new += 1e-5   
    
        P0 = Pressure[-1]
        m0 = m[-1]
    R = r_array[-1]
    M = m[-1]
    compactness = R/M
    return R,M,r_array,P_array,m_array,compactness

