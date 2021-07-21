from TOV import TOV
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functions import rel, mass, y0
from EoS import EoSclass
from functions import tauf1, minP
from interPwEoS import rhof, rhofint

def AllTOV(P0, rho1, rho2=0, r_c=np.Inf, y1=0, y2=0, order=["incom", "incom"]):
  
        #r_c = radius of transition to second fluid
        #y = interaction strength 
        #order: list of types of fluid to be used for fluid 1 and 2.
        #use "incom" for incompressible, "nfermi" for non-interacting, and "intfermi" for interacting

    def TOV(P0,rho_0,tau, r_c, r_initial,m_initial):

        def diff(x,r):
            rho = rho_0
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
        
        x0 = [P0,m_initial]
        r_new = r_initial
        t_span =np.linspace(r_new,tau+r_new,10)

        
        P_array  = np.array([P0])
        r_array  = np.array([r_new])
        m_array  = np.array([0])

        while True:
            sol = odeint(diff,x0,t_span)
            P = sol[:,0] 
            m = sol[:,1]

            P_array = np.append(P_array,P)
            print(P_array[-1])
            r_array = np.append(r_array,t_span)
            m_array = np.append(m_array,m)
            
            if (P <= limit).any():
                index = np.where(P<= limit)
                i = index[0][0]
                if i ==0:
                    R = r_array[-1]
                    M = m_array[-1]
                else:
                    R = t_span[i-1]
                    M = m[i-1]
                compactness = R/M
                print('Star found with R= ',R,'& M=',M, 'Compactness(R/M > 2.25) = ',(R)/(M),'(1 Step Profile) with',len(r_array),'steps')
                break
           
            elif r_array[-1] > r_c and first:
                R = r_array[-1]
                M = m_array[-1]
                print(R)
                print(M)
                compactness = R/M
                break

                    
            t_span = np.linspace(t_span[-1],tau+t_span[-1],10)
            x0 = [P[-1],m[-1]]

        R1 = R
        M1 = M
        comp = compactness
            
        return R1,M1,r_array,P_array,m_array,comp
    

    E = EoSclass

    def TOVEoS(P0, r_c, r_initial, m_initial):

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
        x0 = [P0,m_initial]
        int_P = P0
        limit = minP(int_P)

        tau = tauf1(P0)
        r_new = r_initial
        t_span =np.linspace(r_new,tau+r_new,10)
        
        P_array  = np.array([P0])
        r_array  = np.array([r_new])
        m_array  = np.array([m_initial])

        while True:
            # print(rho)
        
            sol = odeint(diff,x0,t_span)
            P = sol[:,0] 
            m = sol[:,1]

            print(P, "fermi")
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
            elif r_array[-1] > r_c and first:
                R = r_array[-1]
                M = m_array[-1]
                print(R)
                print(M)
                compactness = R/M
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

    def TOVEoSInt(P0,y,r_c,r_initial,m_initial):

    # E = EoSIntClass(y)

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
        
        rho = rhofint(P0,y)
        x0 = [P0,m_initial]
        int_P = P0
        limit = minP(int_P)

        tau = tauf1(P0)
        r_new = r_initial
        t_span =np.linspace(r_new,tau+r_new,10)
        
        P_array  = np.array([])
        r_array  = np.array([])
        m_array  = np.array([])

        while True:
            # print(rho)
        
            sol = odeint(diff,x0,t_span)
            P = sol[:,0] 
            m = sol[:,1]

            print(P[-1],rho)
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
            rho = rhofint(P[-1],y)
            tau = tauf1(P[-1])
            t_span = np.linspace(t_span[-1],tau+t_span[-1],10)
            x0 = [P[-1],m[-1]]

        R1 = R
        M1 = M
        comp = compactness
            
        return R1,M1,r_array,P_array,m_array,comp


    rho = [rho1,rho2]
    first = True
    limit = P0*(1e-6)

    

    if order[0] == "incom":
        sol = TOV(P0, rho[0], 1e-6, r_c, 1e-10, 0)
        r1 = sol[2]
        P1 = sol[3]
        m1 = sol[4]
    elif order[0] == "nfermi":
        sol = TOVEoS(P0, r_c, 1e-10, 0)
        r1 = sol[2]
        P1 = sol[3]
        m1 = sol[4]
    elif order[0] == "intfermi":
        sol = TOVEoSInt(P0, y1, r_c, 1e-10, 0)
        r1 = sol[2]
        P1 = sol[3]
        m1 = sol[4]
    else:
        print("please enter a valid command")
    
    first = False
    if order[1] == "incom":
        sol = TOV(P1[-1], rho[1], 1e-6, r_c, r1[-1], m1[-1])
        r2 = sol[2]
        P2 = sol[3]
        m2 = sol[4]
    elif order[1] == "nfermi":
        sol = TOVEoS(P1[-1], r_c, r1[-1], m1[-1])
        r2 = sol[2]
        P2 = sol[3]
        m2 = sol[4]
    elif order[1] == "intfermi":
        sol = TOVEoSInt(P1[-1], y2, r_c, r1[-1], m1[-1])
        r2 = sol[2]
        P2 = sol[3]
        m2 = sol[4]
    r = np.append(r1, r2)
    P = np.append(P1, P2)
    m = np.append(m1, m2)

    return r, P, m, r1, r2, P1, P2, m1, m2