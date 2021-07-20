import numpy as np
from scipy.integrate import odeint
from EoS import EoSclass
from functions import findXPoint,tauf1
from scipy.interpolate import interp1d

E = EoSclass

def rhof(P):
    i = min(E.EoSP,key=lambda x:abs(x-P)) #find closest value to P0
    a = np.where(E.EoSP==i) # index of closest point to P0
    index =a[0][0] #index of closest P0 (a outputs 2 dim. array)
    f = interp1d([E.EoSP[index-1],E.EoSP[index],E.EoSP[index+1]],[E.EoSrho[index-1],E.EoSrho[index],E.EoSrho[index+1]])    
    return f(P)

def TOVEoS(P0):

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
    int_P = P0
    limit = int_P*(1e-3)

    tau = tauf1(P0)
    r_new = 1e-10
    t_span =np.linspace(r_new,tau+r_new,10)
    
    P_array  = np.array([])
    r_array  = np.array([])
    m_array  = np.array([])

    while True:
        # print(rho)
    
        sol = odeint(diff,x0,t_span)
        P = sol[:,0] 
        m = sol[:,1]

        # print(P)
        if (P <= limit).any():
            index = np.where(P<= limit)
            i = index[0][0]
            if i ==0:
                R = r_array[-1][-1]
                M = m_array[-1][-1]
            else:
                R = t_span[i-1]
                M = m[i-1]
            compactness = R/M
            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M) = ',(R)/(M),'(1 Step Profile) with',len(r_array),'steps')
            break

        P_array = np.append(P_array,P)
        r_array = np.append(r_array,t_span)
        m_array = np.append(m_array,m)
        rho = rhof(P[-1])
        tau = tauf1(P[-1])
        t_span = np.linspace(t_span[-1],tau+t_span[-1],10)
        x0 = [P[-1],m[-1]]

    R1 = R
    M1 = M
    comp = compactness
        
    return R1,M1,r_array,P_array,m_array,comp

# TOV test
# import matplotlib.pyplot as plt
# sol = TOVEoS(1e-7)
# plt.plot(sol[2],sol[3])
# plt.show()