import numpy as np
from functions import EosP, EosRho,intTerm
from scipy.integrate import quad

class EoSclass():
    
    eps = np.finfo(float).eps
    
    kmin = 0
    # kmin = eps
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
        EoSP1 = np.append(EoSP1,a[0])
        EoSrho1 = np.append(EoSrho1,b[0])

    EoSP = np.unique(EoSP1)
    EoSrho = np.unique(EoSrho1)
# import matplotlib.pyplot as plt
# plt.loglog(EoSclass.EoSrho,EoSclass.EoSP,color='red')
# print('min=', min(EoSclass.EoSP),'max=',max(EoSclass.EoSP))
# print('min=', min(EoSclass.EoSrho),'max=',max(EoSclass.EoSrho))
# plt.show()


class EoSIntClass():

    def __init__(self,y):
        self.y = y
    
        eps = np.finfo(float).eps
    
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
            EoSP1 = np.append(EoSP1,a[0]+intTerm(self.y,kmax[i]))
            EoSrho1 = np.append(EoSrho1,b[0]+intTerm(self.y,kmax[i]))

        self.P = np.unique(EoSP1)
        self.rho = np.unique(EoSrho1)



# np.savetxt("eosdata.txt", kmax, delimiter=",")