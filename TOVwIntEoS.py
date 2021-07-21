import numpy as np
from scipy.integrate import quad
from scipy.integrate import odeint
from functions import tauf1,minP,findXPoint,tauR

pi = np.pi
def EosRho(k):
    m_f =1
    return (1/(pi**2))*(k**2)*(np.sqrt(m_f**2 + k**2))
    
def EosP(k):
    m_f =1
    return (1/(3*pi**2))*((k**4)/(np.sqrt(m_f**2 + k**2)))

def intTerm(y,z):
    return (((1/(3*pi**2)))**2)*(y**2)*(z**6)

def TOVintEoS(P0,y):

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
        b =quad(EosRho,kmin,kmax[i]) 
        EoSP1 = np.append(EoSP1,a[0]+intTerm(y,kmax[i]))
        EoSrho1 = np.append(EoSrho1,b[0]+intTerm(y,kmax[i]))

    EoSP = np.unique(EoSP1)
    EoSrho = np.unique(EoSrho1)

    def rhof(P_v):            # Pressure on x-axis and Rho on y-axis
        i = min(EoSP,key=lambda x:abs(x-P_v)) #find closest value to P0
        a = np.where(EoSP==i) # index of closest point to P0
        index =a[0][0] #index of closest P0 (a outputs 2 dim. array)
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
    limit = minP(int_P)

    # tau = tauR(rho)
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

        # print(P[-1],rho)

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
        # tau = tauR(rho)
        tau = tauf1(P[-1])
        t_span = np.linspace(t_span[-1],tau+t_span[-1],10)
        x0 = [P[-1],m[-1]]

    R1 = R
    M1 = M
    comp = compactness
        
    return R1,M1,r_array,P_array,m_array,comp

# TOV test
# import matplotlib.pyplot as plt
# sol = TOVintEoS(10,1000)
# plt.plot(sol[2],sol[3])
# plt.show()
