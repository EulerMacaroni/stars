import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from functions1 import mass, Mass,Nrel,convED, rel,convM #, dP_dr
# define variables
eps = np.finfo(float).eps
# define variables
def TOV(P0,rho):
    eps = np.finfo(float).eps
    rho_0 = rho 
    P = P0
    def dP_dr(P,r):
        k  = mass_step * rho_0/(r**2)
        k1 = 1 + P/rho_0
        k2 = 1 + (4*np.pi * r**3 * P )/ (mass_step)
        k3 = 1 / (1 - (2*mass_step) / (r))
        dPdr = -(k)*(k1*k2*k3)
        return dPdr
    def dm_dr(M,r):
        dmdr = 4*np.pi * r**2 * rho_0
        return dmdr
    # # solving Mass and Pressure differential equations
    P_array = np.array([])
    r_array = np.array([])
    m_array = np.array([])
    r_new = eps
    # mass_step = 0
    m0 = eps
    # print(P)
    r_step_sz = 0.001
    while True:
        # mass_step = mass_step + dm_dr(r_new)*0.01
        mass = solve_ivp(dm_dr,(r_new,r_new,r_step_sz),[m0])
        mass_step = mass.y[0][-1]
        P = solve_ivp(dP_dr, (r_new,r_new+r_step_sz), [P])
        # P_array = np.append(P_array, P.y[0][-1])
        P_array = np.append(P_array, P.y[0])
        r_array  = np.append(r_array, P.t)
        m_array = np.append(m_array,mass.y[0][-1])
        if P.y[0][-1] <= 1e-10:
            R = r_array[-1]
            M = mass_step
            break
        r_new += r_step_sz
        P = P.y[0][-1]
        m0 = mass.y[0][-1]
    R = r_array[-1]
    M = mass_step
    compactness = M/R
    return R, M, compactness, P_array,r_array,m_array


def two_fluid_TOV(P0, rho_0, rho_1=0):
    if not rho_1:
        rho_1 = rho_0/2
    
    P = P0
    
    def dP_dr(P,r,rho):
        k  = mass_step * rho/(r**2)
        k1 = 1 + P/rho
        k2 = 1 + 4*np.pi * r**3 * P / mass_step
        k3 = 1 / (1 - 2*mass_step / r)
        dPdr = -(k)*(k1*k2*k3)
        return dPdr
    def dm_dr(r, rho_0):
        dmdr = 4*np.pi * r**2 * rho_0
        return dmdr
    # # solving Mass and Pressure differential equations for the first fluid
    P_array = np.array([])
    r_array = np.array([])
    m_array = np.array([])
    r_array_norm = np.array([]) # same length as the mass array
    r_new = eps
    mass_step = 0
    r_step_sz = 0.001
    second_fluid = False
    while True:
        mass_step = mass_step + dm_dr(r_new, rho_0)*r_step_sz
        #mass = solve_ivp()
        P = solve_ivp(dP_dr, (r_new,r_new+r_step_sz), [P], args=(rho_0,))
        P_array = np.append(P_array, P.y[0])
        r_array  = np.append(r_array, P.t)
        r_array_norm = np.append(r_array_norm, P.t[-1])
        m_array = np.append(m_array, mass_step)
        if P.y[0][-1] <= 1e-10:
            R = r_array[-1]
            M = mass_step
            P = P.y[0][-1]
            break
        r_new += r_step_sz
        P = P.y[0][-1]

    
    P_array = P_array[:int(len(r_array)/2)]
    r_array = r_array[:int(len(r_array)/2)]
    r_array_norm = r_array_norm[:int(len(m_array)/2)]
    mass_step = m_array[int((len(r_array))/2)]
    m_array = m_array[:int(len(r_array)/2)]

    

    P = P_array[-1]
    r_new = r_array[-1] + r_step_sz

    P_array1 = np.array([])
    r_array1 = np.array([])
    m_array1 = np.array([])

    #solving for the second fluid
    while True:
        mass_step = mass_step + dm_dr(r_new, rho_1)*r_step_sz
        P = solve_ivp(dP_dr, (r_new,r_new+r_step_sz), [P], args= (rho_1,))
        P_array1 = np.append(P_array1, P.y[0])
        r_array1 = np.append(r_array1, P.t)
        r_array_norm = np.append(r_array_norm, P.t[-1])
        m_array1 = np.append(m_array1, mass_step)
        if P.y[0][-1] <= 1e-10:
            R = r_array[-1]
            M = mass_step
            break
        r_new += r_step_sz
        P = P.y[0][-1]


    
    r_array_total = np.append(r_array, r_array1)
    P_array_total = np.append(P_array, P_array1)
    m_array_total = np.append(m_array, m_array1)


    R = r_array[-1]
    M = mass_step
    compactness = M/R
    #returns [0] total radius, [1] total mass, [2] compactness, [3] List of Pressures for both fluids
    # [4] radius list for both fluids, [5] mass list for both fluids, [6] Pressure list for rho_0
    # [7] radius list for rho_0, [8] mass list for rho_0, [9] Pressure list for rho_1, [10] radius list for rho_1
    # [11] mass list for rho_1, [12] for radius list of equal length to the mass list
    return R, M, compactness, P_array_total , r_array_total, m_array_total, P_array, r_array, m_array, P_array1, r_array1, m_array1, r_array_norm




#P vs M/R Plot for one fluid
# rho = 10e-25
# P = np.arange(eps,1,0.001)
# print(len(P))
# r_m = []
# for i in range(len(P)):
#    sol = TOV(P[i],rho)
#    r_m.append(sol[2])
# plt.plot(r_m,P)
# #plt.axvline(x=9/4)
# plt.show()

#P vs M/R Plot for two fluids
rho = 10e-25
P = np.arange(eps,1,0.01)
print(len(P))
r_m = np.array([])
for i in range(len(P)):
   sol = TOV(P[i],rho)
   r_m = np.append(r_m, sol[2])
plt.plot(r_m,P)
plt.ylabel("Pressure $P$")
plt.xlabel("Radius / Mass")
plt.show()

# One Fluid P vs r plot
# sol = TOV(1,10e-25)
# R = sol[0]
# M = sol[1]
# P = sol[3]
# r = sol[4]
# # print("compactness: ", compactness)
# # print("Total Mass: ", M)
# # print("radius: ", R)
# plt.yscale('log')
# plt.plot(r, P)
# plt.legend(["calculated"])
# plt.xlabel('Radius $r$')
# plt.ylabel('Pressure $P$')
# plt.title("One Fluid Solution")
# plt.show()


#two fluid (Pressure and Mass vs. radius)
rho_0 = 10e-25
sol = two_fluid_TOV(1, rho_0, rho_0/2)
R = sol[0]
M = sol[1]
r = sol[12]
m = sol[5]
P0 = sol[6]
r0 = sol[7]
P1 = sol[9]
r1 = sol[10]
# print("compactness: ", compactness)
# print("Total Mass: ", M)
# print("radius: ", R)
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_yscale('log')
ax1.plot(r0, P0)
ax1.plot(r1, P1)
ax2.plot(r, m, color="red")
ax1.legend(["rho_0", "rho_1"], loc="center left")
ax2.legend(["mass"], loc="upper left" )
plt.xlabel('Radius $r$')
ax1.set_ylabel('Pressure $P$')
ax2.set_ylabel("Mass $M$")
plt.title("Two Fluid Solution")
plt.show()