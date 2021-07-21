import re
import numpy as np
from EoS import EoSclass,EoSIntClass
from functions import findXPoint
from scipy.interpolate import interp1d

E = EoSclass
def rhof(P):            # Pressure on x-axis and Rho on y-axis
    i = min(E.EoSP,key=lambda x:abs(x-P)) #find closest value to P0
    a = np.where(E.EoSP==i) # index of closest point to P0
    index =a[0][0] #index of closest P0 (a outputs 2 dim. array)
    p1 = E.EoSP[index]
    rho1 = E.EoSrho[index]

    if P > E.EoSP[index]:
        p2 = E.EoSP[index+1]
        rho2 = E.EoSrho[index+1]
        f1 = findXPoint(rho1,rho2,p1,p2,P)
        # f2 = interp1d([p1,p2],[rho1,rho2])

    elif P < E.EoSP[index]:
        p2 = E.EoSP[index-1]
        rho2 = E.EoSrho[index-1]
        f1 = findXPoint(rho2,rho1,p2,p1,P)
        # f2 = interp1d([p2,p1],[rho2,rho1])
    elif P == E.EoSP[index]:
        f1 = E.EoSrho[index]

    return f1


def Pforrho(rho):
    i = min(E.EoSrho,key=lambda x:abs(x-rho)) #find closest value to P0
    a = np.where(E.EoSrho==i) # index of closest point to P0
    index =a[0][0] #index of closest P0 (a outputs 2 dim. array)

    p1 = E.EoSP[index]
    rho1 = E.EoSrho[index]

    if rho > E.EoSrho[index]:
        p2 = E.EoSP[index+1]
        rho2 = E.EoSrho[index+1]
        f = findXPoint(p1,p2,rho1,rho2,rho)
        # f2 = interp1d([p1,p2],[rho1,rho2])

    elif rho < E.EoSrho[index]:
        p2 = E.EoSP[index-1]
        rho2 = E.EoSrho[index-1]
        f = findXPoint(p2,p1,rho2,rho1,rho)

    elif rho == E.EoSrho[index]:
        f = E.EoSP[index]

    return f