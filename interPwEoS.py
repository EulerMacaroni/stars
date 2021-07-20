import numpy as np
from EoS import EoSclass
from scipy.interpolate import interp1d

E = EoSclass

def rhof(P):            # Pressure on x-axis and Rho on y-axis
    i = min(E.EoSP,key=lambda x:abs(x-P)) #find closest value to P0
    a = np.where(E.EoSP==i) # index of closest point to P0
    index =a[0][0] #index of closest P0 (a outputs 2 dim. array)

    if P > E.EoSP[index]:
        p1 = E.EoSP[index]
        rho1 = E.EoSrho[index]
        p2 = E.EoSP[index+1]
        rho2 = E.EoSrho[index+1]
        f = interp1d([p1,p2],[rho1,rho2])

    else:
        p1 = E.EoSP[index]
        rho1 = E.EoSrho[index]
        p2 = E.EoSP[index-1]
        rho2 = E.EoSrho[index-1]
        f = interp1d([p2,p1],[rho2,rho1])

    return f(P)