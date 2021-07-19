import numpy as np
from functions import EosP, EosRho,intTerm
from scipy.integrate import quad

class EoSclass():
    
    eps = np.finfo(float).eps
    
    kmin = eps
    # kmin = eps
    kmax = np.linspace(kmin+eps,9000,10000)

    EoSrho = np.array([])
    EoSP = np.array([])
    for i in range(len(kmax)):
        a =quad(EosP,kmin,kmax[i]) 
        b =quad(EosRho,kmin,kmax[i]) 
        EoSP = np.append(EoSP,a[0])
        EoSrho = np.append(EoSrho,b[0])

# import matplotlib.pyplot as plt
# plt.plot(EoSclass.EoSP,EoSclass.EoSrho,'-.')
# print('min=', min(EoSclass.EoSP),'max=',max(EoSclass.EoSP))
# print('min=', min(EoSclass.EoSrho),'max=',max(EoSclass.EoSrho))
# plt.show()


class EoSIntClass():

    def __init__(self,y):
        self.y = y
    
        eps = np.finfo(float).eps
    
        kmin = eps
        # kmin = eps
        kmax = np.linspace(kmin+eps,9000,10000)

        EoSrho = np.array([])
        EoSP = np.array([])
        for i in range(len(kmax)):
            a =quad(EosP,kmin,kmax[i]) 
            b =quad(EosRho,kmin,kmax[i]) 
            EoSP = np.append(EoSP,a[0]+intTerm(self.y,kmax[i]))
            EoSrho = np.append(EoSrho,b[0]+intTerm(self.y,kmax[i]))
        self.P = EoSP
        self.rho = EoSrho