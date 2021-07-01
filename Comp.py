import numpy as np
import matplotlib.pyplot as plt
from TOV import TOV
from functions import rel, mass

# eps = np.finfo(float).eps
# P = np.linspace(1e-5,1e3,1000)
# rho = 261            # (107-108)
# r_m = np.array([])
# for i in range(len(P)):
#     sol = TOV(P[i],rho)
#     r_m = np.append(r_m,sol[5])

# plt.ylabel('Central Pressure, $P_0$')
# plt.xlabel('R/M')
# plt.xlim([0,20])
# plt.title('Central Pressure vs R/M with rho = %f'%rho)
# plt.axvline(x=9/4,label='Buchdahl limit',color='green',linestyle = 'dashed')
# plt.plot(r_m,P)
# plt.show()

# h = 1000
# # rho = 1e4
# while True:
#     P = 1e6
#     rho = 1000+h
#     print('P0 = ',P,'and rho=',rho)
#     sol = TOV(P,rho)
#     if sol[5] <= 2.25:
#         print('illegal star')
#         # break
#     h += 1000