import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint
from functions import EosP, EosRho
from scipy.interpolate import interp1d
from scipy.integrate import quad

class EoSclass():
    
    eps = np.finfo(float).eps
    
    kmin = 1e-8
    kmax = np.linspace(kmin+eps,9000,90000)

    EoSrho = np.array([])
    EoSP = np.array([])
    for i in range(len(kmax)):
        a =quad(EosP,kmin,kmax[i]) 
        b =quad(EosRho,kmin,kmax[i]) 
        EoSP = np.append(EoSP,a[0])
        EoSrho = np.append(EoSrho,b[0])