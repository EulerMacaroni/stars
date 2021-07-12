import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions import rel, mass,y0
from TOV import TOV

# single density profile TOV solver

def TOVstep(m_r,rho,tau):

    sol = TOV(m_r,rho,tau)
    P = y0(m_r,rho) 
    r_c = sol[0]/2

    def rho_0(r):
        if r <= r_c:
            return rho
        else:
            return rho/2

    def dm_dr(r,M):
        return (4)*(np.pi)*(r**2)*rho_0(r)

    def dP_dr(r,P):
        k  = (rho_0(r))/(r)
        k1 = 1 + (P/rho_0(r))
        k2 = mstep + ((4*np.pi*(r**3)*P))
        k3 = (r - (2*mstep))**(-1)
        dPdr = -(k)*(k1*k2*k3)
        return dPdr

    eps = np.finfo(float).eps
    P0 = P
    r_new = 1e-10
    m0 = 0
    t_span =np.linspace(r_new,tau+r_new,50)
    
    # # solving Mass and Pressure differential equations

    P_array  = np.array([])
    r_array  = np.array([])
    m_array  = np.array([])

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
            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M > 2.25) = ', (R)/(M),'(2 Step Profile) with', len(r_array),'steps')
            break
        
        t_span = np.linspace(r_array[-1],tau+r_array[-1],50)
    
        P0 = Pressure[-1]
        m0 = m[-1]
    R = r_array[-1]
    M = m[-1]
    compactness = R/M
    return R,M,r_array,P_array,m_array,compactness, r_c,P