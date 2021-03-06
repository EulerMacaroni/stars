import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions import rel, mass, y0
from TOV import TOV


def TOV2step(P0,rho_0,tau):
    
    sol = TOV(P0,rho_0,tau) 
    r_c = sol[0]/2

    def rho(r):
        if r <= r_c:
            return rho_0
        else:
            return rho_0/2

    def diff(x,r):
        P = x[0]
        m = x[1]
        # TOV
        k  = (rho(r))/(r)
        k1 = 1 + (P/rho(r))
        k2 = m + ((4*np.pi*(r**3)*P))
        k3 = (r - (2*m))**(-1)
        dPdr = -(k)*(k1*k2*k3) 
        # Mass
        dmdr = 4*np.pi*(r**2)*rho(r)
        return [dPdr,dmdr]

    eps = np.finfo(float).eps
    
    x0 = [P0,0]
    
    r_new = 1e-10
    t_span =np.linspace(r_new,tau+r_new,10)
    
    P_array  = np.array([])
    r_array  = np.array([])
    m_array  = np.array([])

    while True:
    
        sol = odeint(diff,x0,t_span)
        P = sol[:,0] 
        m = sol[:,1]
        
        if (P <= 1e-10).any():
            index = np.where(P<= 1e-10)
            i = index[0][0]
            if i ==0:
                R = r_array[-1][-1]
                M = m_array[-1][-1]
            else:
                R = t_span[i-1]
                M = m[i-1]
            compactness = R/M
            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M > 2.25) = ', (R)/(M),'(2 Step Profile) with', len(r_array),'steps')
            break

        P_array = np.append(P_array,P)
        r_array = np.append(r_array,t_span)
        m_array = np.append(m_array,m)
                
        t_span = np.linspace(t_span[-1],tau+t_span[-1],10)
        x0 = [P[-1],m[-1]]

    R1 = R
    M1 = M
    comp = compactness
        
    return R1,M1,r_array,P_array,m_array,comp,r_c
