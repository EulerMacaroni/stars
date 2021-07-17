from warnings import resetwarnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from EoS import EoSclass
from functions import rel, mass, y0,findXPoint


def TOVEoS(P0,tau):

    E = EoSclass

    def rhof(P):
        i = min(E.EoSP,key=lambda x:abs(x-P)) #find closest value to P0
        a = np.where(E.EoSP==i) # index of closest point to P0
        index =a[0][0] #index of closest P0 (a outputs 2 dim. array)
    
        x2 = E.EoSrho[index+1]
        x1 = E.EoSrho[index]
        y2 = E.EoSP[index+1]
        y1 = E.EoSP[index]

        x3 = findXPoint(x1,x2,y1,y2,P)
        return x3


    def diff(x,r):
        P = x[0]
        m = x[1]
        # TOV
        k  = (rho)/(r)
        k1 = 1 + (P/rho)
        k2 = m + ((4*np.pi*(r**3)*P))
        k3 = (r - (2*m))**(-1)
        dPdr = -(k)*(k1*k2*k3) 
        # Mass
        dmdr = 4*np.pi*(r**2)*rho
        return [dPdr,dmdr]

    eps = np.finfo(float).eps
    
    rho = rhof(P0)
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
            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M > 2.25) = ',(R)/(M),'(1 Step Profile) with',len(r_array),'steps')
            break

        P_array = np.append(P_array,P)
        r_array = np.append(r_array,t_span)
        m_array = np.append(m_array,m)
                
        t_span = np.linspace(t_span[-1],tau+t_span[-1],10)
        x0 = [P[-1],m[-1]]
        rho = rhof(P[-1])

    R1 = R
    M1 = M
    comp = compactness
        
    return R1,M1,r_array,P_array,m_array,comp
