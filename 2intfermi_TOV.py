import numpy as np
from numpy.core.function_base import linspace
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import quad
from functions import findXPoint, minP
import numpy as np

# we want core made of dark matter (interacting fermi EoS) and rest of the star is made of quark matter or fermi (diff interaction)
# define the core by density 

pi = np.pi
B = (145)**4
m_f = 1
Mp = 1.3394 * (10**-1)
# G = Mp**(-2)
G =1
a = (Mp**3)/(m_f**4)
b = (Mp)/(m_f**2)
x1 = np.sqrt(B)
x2 = 4*B
x3 = m_f**4

def EosRho(k):
    return (1/(pi**2))*(k**2)*(np.sqrt(m_f**2 + k**2))
    
def EosP(k):
    return (1/(3*pi**2))*((k**4)/(np.sqrt(m_f**2 + k**2)))

def intTerm(y,z):
    return (((1/(3*pi**2)))**2)*(y**2)*(z**6)

def EoS2(rho_0, y1 ,y2):
    kmin = 0
    kmax = np.array([])
    k1 = np.linspace(1e-16,1e-15,500,endpoint=True)

    for i in range(0,21):
        k = k1*(10**i)
        kmax = np.append(kmax,k)

    EoSrho_1 = np.array([])
    EoSP_1 = np.array([])
    EoSrho_2 = np.array([])
    EoSP_2 = np.array([])

    for i in range(len(kmax)):
        a =quad(EosP,kmin,kmax[i]) 
        b =quad(EosRho,kmin,kmax[i])            # if u cut at 1e-2, lin (1e-2 ,1000) --> P 1/3 
        EoSP_1 = np.append(EoSP_1,a[0]+intTerm(y1,kmax[i]))
        EoSrho_1 = np.append(EoSrho_1,b[0]+intTerm(y1,kmax[i]))
        EoSP_2 = np.append(EoSP_2,a[0]+intTerm(y2,kmax[i]))
        EoSrho_2 = np.append(EoSrho_2,b[0]+intTerm(y2,kmax[i]))

    EoSP1 = np.unique(EoSP_1)
    EoSrho1 = np.unique(EoSrho_1)
    EoSP2 = np.unique(EoSP_1)
    EoSrho2 = np.unique(EoSrho_1)

    def Rho4P(P_v,a):            
        if a ==1:
            EoSP = EoSP1
            EoSrho = EoSrho1
        elif a==2:
            EoSP = EoSP2
            EoSrho = EoSrho2

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

    def P4Rho(rho,b):
        if b ==1:
            EoSP = EoSP1
            EoSrho = EoSrho1
        elif b==2:
            EoSP = EoSP2
            EoSrho = EoSrho2

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

    def diff(x,r):
        P = x[0]
        m = x[1]

        k  = (rho)/(r)
        k1 = 1 + (P/rho)
        k2 = m + ((4*pi*(r**3)*P))
        k3 = (r - (2*m))**(-1)

        dpdr = -(k)*(k1*k2*k3) 
        dmdr = 4*np.pi*(r**2)*rho
        return [dpdr,dmdr]

    def choose(rho):
        if rho >= 1e-2:
            return 1
        else:
            return 2

    rho = rho_0
    P0 = P4Rho(rho,choose(rho))

    x0 = [P0,0]
    int_P = P0
    limit = minP(int_P)

    tau = 1e-10
    r_new = 1e-10
    t_span =np.linspace(r_new,tau+r_new,10)
    
    P_array  = np.array([])
    r_array  = np.array([])
    m_array  = np.array([])
    rho_array = np.array([])

    while True:
        # print(rho)
    
        sol = odeint(diff,x0,t_span)
        P = sol[:,0] 
        m = sol[:,1]

        P_array = np.append(P_array,P)
        r_array = np.append(r_array,t_span)
        m_array = np.append(m_array,m)
        rho_array = np.append(rho_array,rho)

        # print(P)
        if (P <= limit).any():
            index = np.where(P_array<= limit)
            i = index[0][0]
            R = r_array[i]
            M = m_array[i]
            P_val = P_array[0:i]
            r_val = r_array[0:i]
            m_val = m_array[0:i]
            rho_val = rho_array[0:i]

            comp = R/M
            print('Star found with R= ',R,'& M=',M, 'Compactness(R/M) = ',(R)/(M),'(2 Step Profile) with',len(r_array),'steps')
            break

        rho = Rho4P(P[-1],choose(rho_array[-1]))
        tau = tau*1.01
        t_span = np.linspace(t_span[-1],tau+t_span[-1],10)
        x0 = [P[-1],m[-1]]

        
    return R,M,r_val,P_val,m_val,comp

y1 = 10
y2 = 1e-2

rho1 = np.linspace(10**(-8),10**(-7),10,endpoint=True)
rho = np.array([])
for i in range(0,9):
    k = rho1*(10**i)
    rho = np.append(rho,k)

R = np.array([])
M = np.array([])
comp = np.array([])

for i in range(len(rho)):
    sol = EoS2(rho[i], y1 ,y2)
    R = np.append(R,sol[0])
    M = np.append(M,sol[1])
    comp = np.append(comp,sol[5])

fig1 = plt.figure(1)
plt.ylabel('dimensionaless M')
plt.xlabel('dimensionaless R')
plt.plot(R, M,'.' ,color='black')

fig1 = plt.figure(2)
plt.plot(R, np.power(comp,-1),'.')
plt.axhline(y=4/9,color='green',linestyle = 'dashed')
plt.legend(['2 Step','Buchdahl limit'])
plt.xlabel('R')
plt.ylabel('M/R')
plt.ylim([0,4/9 +0.1 ])
plt.grid()

plt.show()
