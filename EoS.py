import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint
from functions import EosP, EosRho
from scipy.interpolate import interp1d
from scipy.integrate import quad

class EoS():
    
    eps = np.finfo(float).eps
    
    kmin = 1e-8
    kmax = np.linspace(kmin+eps,5000,10000)

    EoSrho = np.array([])
    EoSP = np.array([])
    for i in range(len(kmax)):
        a =quad(EosP,kmin,kmax[i]) 
        b =quad(EosRho,kmin,kmax[i]) 
        EoSP = np.append(EoSP,a[0])
        EoSrho = np.append(EoSrho,b[0])



# Plot of pressure vs rho 

plt.loglog(EoS.EoSrho,EoS.EoSP)
plt.xlabel('$\epsilon$')
plt.ylabek('$P$')
plt.title('pressure vs density')
plt.legend(['Pressure'])
plt.show()

