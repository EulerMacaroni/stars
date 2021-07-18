import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint
from functions import EosP, EosRho,intTerm
from scipy.interpolate import interp1d
from scipy.integrate import quad

class EoSclass():
    
    eps = np.finfo(float).eps
    
    kmin = 1e-8
    # kmin = eps
    kmax = np.linspace(kmin+eps,5000,10000)

    EoSrho = np.array([])
    EoSP = np.array([])
    for i in range(len(kmax)):
        a =quad(EosP,kmin,kmax[i]) 
        b =quad(EosRho,kmin,kmax[i]) 
        EoSP = np.append(EoSP,a[0])
        EoSrho = np.append(EoSrho,b[0])

# print('min=', min(EoSclass.EoSP),'max=',max(EoSclass.EoSP))

class EoSIntClass():

    def __init__(self,y):
        self.y = y
    
        eps = np.finfo(float).eps
    
        kmin = 1e-8
        # kmin = eps
        kmax = np.linspace(kmin+eps,5000,10000)

        EoSrho = np.array([])
        EoSP = np.array([])
        for i in range(len(kmax)):
            a =quad(EosP,kmin,kmax[i]) 
            b =quad(EosRho,kmin,kmax[i]) 
            EoSP = np.append(EoSP,a[0]+intTerm(self.y,kmax[i]))
            EoSrho = np.append(EoSrho,b[0]+intTerm(self.y,kmax[i]))
        self.P = EoSP
        self.rho = EoSrho