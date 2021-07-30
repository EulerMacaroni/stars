import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import quad
from functions import findXPoint, minP

pi = np.pi
B = (145)**4
m_f = 1
def EosRho(k):
    return (1/(pi**2))*(k**2)*(np.sqrt(m_f**2 + k**2))
    
def EosP(k):
    return (1/(3*pi**2))*((k**4)/(np.sqrt(m_f**2 + k**2)))

def intTerm(y,z):
    return (((1/(3*pi**2)))**2)*(y**2)*(z**6)


def TOVmix(rho_01,rho_02,y):

    kmin = 0
    kmax = np.array([])
    k1 = np.linspace(1e-16,1e-15,500,endpoint=True)

    for i in range(0,21):
        k = k1*(10**i)
        kmax = np.append(kmax,k)

    EoSrho1 = np.array([])
    EoSP1 = np.array([])

    for i in range(len(kmax)):
        a =quad(EosP,kmin,kmax[i]) 
        b =quad(EosRho,kmin,kmax[i])            # if u cut at 1e-2, lin (1e-2 ,1000) --> P 1/3 
        EoSP1 = np.append(EoSP1,a[0]+intTerm(y,kmax[i]))
        EoSrho1 = np.append(EoSrho1,b[0]+intTerm(y,kmax[i]))

    EoSP = np.unique(EoSP1)
    EoSrho = np.unique(EoSrho1)

    def Rho4P(P_v):            
        i = min(EoSP,key=lambda x:abs(x-P_v)) 
        a = np.where(EoSP==i) 
        index =a[0][0] 
        p1 = EoSP[index]
        rho1 = EoSrho[index]

        if P_v > EoSP[index]:
            p2 = EoSP[index+1]
            rho2 = EoSrho[index+1]
            f1 = findXPoint(rho1,rho2,p1,p2,P_v)
            # f2 = interp1d([p1,p2],[rho1,rho2])

        elif P_v < EoSP[index]:
            p2 = EoSP[index-1]
            rho2 = EoSrho[index-1]
            f1 = findXPoint(rho2,rho1,p2,p1,P_v)
            # f2 = interp1d([p2,p1],[rho2,rho1])
        elif P_v == EoSP[index]:
            f1 = EoSrho[index]

        return f1

    def P4Rho(rho):
        i = min(EoSrho,key=lambda x:abs(x-rho)) #find closest value to P0
        a = np.where(EoSrho==i) # index of closest point to P0
        index =a[0][0] #index of closest P0 (a outputs 2 dim. array)

        p1 = EoSP[index]
        rho1 = EoSrho[index]

        if rho > EoSrho[index]:
            p2 = EoSP[index+1]
            rho2 = EoSrho[index+1]
            f = findXPoint(p1,p2,rho1,rho2,rho)
            # f2 = interp1d([p1,p2],[rho1,rho2])

        elif rho < EoSrho[index]:
            p2 = EoSP[index-1]
            rho2 = EoSrho[index-1]
            f = findXPoint(p2,p1,rho2,rho1,rho)

        elif rho == EoSrho[index]:
            f = EoSP[index]

        return f

    def mitbagmod4rho(rho):
        return (1/3)*(rho - 3*B)
    
    def mitbagmodel4P(P):
        return 3*P + 4*B


    P01 = mitbagmod4rho(rho_01)
    P02 = P4Rho(rho_02)

    def diff(x,r):
        P1 = x[0]
        P2 = x[1]
        m  = x[2]

        k01 = (rho1)/(r)
        k11 = 1 + (P1/rho1)
        k21 = m + 4*pi*(r**3)*(P1+P2)
        k31 = (r - (2*m))**(-1)

        k02 = (rho2)/(r)
        k12 = 1 + (P2/rho2)
        k22 = m + 4*pi*(r**3)*(P1+P2)
        k32 = (r - (2*m))**(-1)

        dP1dr = -(k01)*(k11*k21*k31) 
        dP2dr = -(k02)*(k12*k22*k32) 
        dmdr  = 4*pi*(r**2)*(rho1+rho2)

        return [dP1dr,dP2dr,dmdr]

    rho1 = rho_01
    rho2 = rho_02

    y0 = [P01,P02,0]

    limit = [minP(P01),minP(P02)]
    tau = 1e-10
    r_int = 1e-10

    P1_array = np.array([])
    P2_array = np.array([])
    r_array  = np.array([])
    m_array  = np.array([])
    rho1_array = np.array([]) 
    rho2_array = np.array([]) 
    

    r = np.linspace(r_int,tau+r_int,20)

    while True:
    
        sol = odeint(diff,y0,r)
        P1  = sol[:,0] 
        P2  = sol[:,1]
        m   = sol[:,2]

        # print('pressure =',P[-1],'density=',rho)

        P1_array = np.append(P1_array,P1)
        P2_array = np.append(P2_array,P2)
        r_array = np.append(r_array,r)
        m_array = np.append(m_array,m)

        if (P1 <= limit[0]).any() and (P2 <= limit[1]).any():

            index1 = np.where(P1_array<= limit[0])
            index2 = np.where(P2_array<= limit[1])

            i1 = index1[0][0]
            i2 = index2[0][0]
            r_vec = np.array([r_array[i1],r_array[i2]])

            if r_vec[0] >= r_vec[1]:  
                Radius = r_vec[0]
                i = i1
            else:
                Radius = r_vec[1]
                i = i2
            
            M = m_array[i]
            R = Radius
            r_val = r_array[0:i]
            m_val = m_array[0:i]
            P1_val = P1_array[0:i]
            P2_val = P2_array[0:i]
            rho1_val = rho1_array[0:i]
            rho2_val = rho2_array[0:i]

            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M) = ',(R)/(M),'(Two EoS) with',len(r_array),'steps')

            break

        rho1 = mitbagmodel4P(P1[-1])
        rho2 = Rho4P(P2[-1])
        tau  = tau*5
        r = np.linspace(r[-1],tau+r[-1],20)
        y0 = [P1[-1],P2[-1],m[-1]]

    return R,M,r_val,m_val,[P1_val,P2_val],[rho1_val,rho2_val]

# sol = TOVmix(1e-3,1e-3,0)
# plt.plot(sol[2],sol[3])
# plt.plot(sol[0],sol[1],'o')
# plt.show()