import numpy as np
import matplotlib.pyplot as plt
from TOV import TOV
from functions import rel, mass
from TOV2step import TOV2step

eps = np.finfo(float).eps
P = np.linspace(1e-5,1e5,100)
rho = 5000           # (107-108)
r_m = np.array([])
for i in range(len(P)):
    print(i)
    sol = TOV(P[i],rho)
    r_m = np.append(r_m,sol[5])
    if sol[5] <= 9/4:
        print('Illegal')

plt.ylabel('Central Pressure, $P_0$')
plt.xlabel('R/M')
plt.xlim([0,20])
plt.title('Central Pressure vs R/M with rho = %f'%rho)
plt.axvline(x=9/4,label='Buchdahl limit',color='green',linestyle = 'dashed')
plt.plot(r_m,P)
plt.show()

# h = 1000
# # rho = 1e4
# while True:
#     P = 1e9
#     rho = 1000+h
#     print('P0 = ',P,'and rho=',rho)
#     sol = TOV(P,rho)
#     if sol[5] <= 2.24:
#         print('illegal star')
#         break
#     h += 1e4



# rho = 1

# # sol1 = TOV2step(1,rho,2.25/2)
# sol = TOV(1,rho)
# plt.yscale('log')
# plt.plot(sol[2],sol[3])
# print('mass calculated=', mass(sol[0],rho))
# print('m_array',sol[4])
# print('sum=',np.sum(sol[4]))
# plt.show()

# rho = 100
# sol = TOV(1e3,rho)
# R = sol[0]
# M = sol[1]
# r_arr = sol[2]

# m_arr = sol[4]
# P = sol[3]
# M2 = mass(R,rho)
# P1 = rel(r_arr,R,M,rho)
# print('analytical mass=',M2, 'and num =',M[0])
# # print('P=',P)
# plt.plot(r_arr,P)
# # print('r=',r_arr)
# plt.plot(r_arr,P1)
# plt.legend(['num','exact'])
# plt.show()


