from interPwEoS import Pforrho
import numpy as np
import matplotlib.pyplot as plt
from TOVwIntEoS import TOVintEoS


# y = [0.01,1,10,10**2,10**3]
y = np.array([])
y1 = np.linspace(0.01,1,5)
for i in range(0,2):
    h = y1*(10**i)
    y = np.append(y,h)
print(y)
rho1 = np.linspace(10**(-8),10**(-7),10,endpoint=True)

rho = np.array([])
for i in range(0,9):
    k = rho1*(10**i)
    rho = np.append(rho,k)

P0 = np.array([])
for i in range(len(rho)):
    P0 = np.append(P0,Pforrho(rho[i]))

# P0 = np.linspace(0.9,4,4,endpoint=True)        # limits = (1e-78,4,10)

R = np.array([])
M = np.array([])
for j in range(len(y)):
    for i in range(len(P0)):
        print('Star #',i,'with P=',P0[i],'and rho=',rho[i],'y=',y[j])
        sol = TOVintEoS(P0[i],y[j])
        # plt.plot(sol[2],sol[4],color='red')
        M = np.append(M,sol[1])
        R = np.append(R,sol[0])

    plt.loglog(R,M,'.',color='black')
# plt.plot(R,M,color='red')

plt.xlim([1,1e4])
plt.ylim([1e-3,200])
plt.xlabel('dimensionaless $R$')
plt.ylabel('dimensionaless $M(R)$')
plt.show()
