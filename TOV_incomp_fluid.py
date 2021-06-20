# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 17:28:11 2021

@author: Benutzer01
"""

import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import pi
#from polytropes import monotrope, polytrope
#from crust import SLyCrust
#from eoslib import get_eos, glue_crust_and_core, eosLib
from scipy.integrate import odeint
#from label_line import label_line

#import palettable as pal
#cmap = pal.colorbrewer.qualitative.Set1_6.mpl_colormap

lss = ['-', '--', '-.', ':']
clr = ['tab:red','tab:blue', 'tab:orange', 'tab:green']

#plot mass
folder = 'data_for_py/'
data1 = []
data1.append(np.genfromtxt(folder + '1outM.dat'))

fig, ax = plt.subplots()

ax.plot(data1[0][:,0], data1[0][:,1], label = "$M_{numerical}$",   linewidth=2,  color=clr[0])
ax.plot(data1[0][:,0], data1[0][:,2], label = "$M_{analytical}$",  linewidth=3, linestyle='dashed',
        color=clr[1])

ax.set_xlabel('$r$',fontsize=12)
ax.set_xlim(0,12)

ax.set_ylabel(r'Mass M',fontsize=12)
ax.set_ylim(0,4.)

ax.set_title(r"The mass of an incompressible fluid")
ax.legend(title="$Mass$", fontsize=12, frameon=False)

plt.minorticks_on()

plt.savefig("plot_M.pdf",dpi=256,transparent=False)

#plot pressure
folder = 'data_for_py/'
data2 = []
data2.append(np.genfromtxt(folder + '1outP.dat'))

fig, bx = plt.subplots()

bx.plot(data2[0][:,0], data2[0][:,1], label = "$P_{numerical}$",   linewidth=2,  color=clr[0])
bx.plot(data2[0][:,0], data2[0][:,2], label = "$P_{analytical}$",  linewidth=3, linestyle='dashed',
        color=clr[1])

bx.set_xlabel('$r$',fontsize=12)
bx.set_xlim(0,12)

bx.set_ylabel(r'Pressure P',fontsize=12)
bx.set_ylim(0,0.004)

bx.set_title(r"The pressure of an incompressible fluid")
bx.legend(title="$Pressure$", fontsize=12, frameon=False)

plt.minorticks_on()

plt.savefig("plot_P.pdf",dpi=256,transparent=False)
# plt.show()
# plt.clf()