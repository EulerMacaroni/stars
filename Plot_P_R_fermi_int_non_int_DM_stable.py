# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 10:35:50 2021

@author: Marie
"""
import matplotlib as mpl  
#mpl.use("pgf")
#mpl.rc('font',family='Times New Roman')
#mpl.rc('font',family='Computer Modern')
#mpl.rc('font',family='Serif')
import matplotlib.pyplot as plt
from matplotlib import rcParams  
from matplotlib import rc 
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

 
lss = ['-', '--', '-.', ':']
clr = ['tab:red','tab:blue', 'tab:orange', 'tab:green' , 'tab:brown' , 'tab:purple']
clr0 = ['css:darkred']
cteal = ['css:teal']
cteal1 = ['xkcd:teal']
cteal2 = ['#528B8B']
clr1 = mcolors.CSS4_COLORS


# change font
#rcParams['mathtext.fontset'] = 'custom'
#rcParams['mathtext.it'] = 'Arial:italic'
#rcParams['mathtext.rm'] = 'Arial'
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif' 
rcParams['font.serif'] = ['Computer Modern']

rcParams['axes.labelsize']=10
#rcParams['axes.labelweight']=600
rcParams['axes.linewidth']=1
rcParams['lines.markersize']=6

rcParams['xtick.labelsize']=14
rcParams['ytick.labelsize']=14

plt.rcParams.update({
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "pgf.preamble": "\n".join([
        # r"\usepackage{url}",            # load additional packages
         r"\usepackage{unicode-math}",   # unicode math setup
         r"\setmainfont{Computer Modern}",  # serif font via preamble
    ])
})

data = []
data.append(np.genfromtxt('1outP0_R_woInt_10m6_wInt_mN_y10_rc1_stable.dat'))
data.append(np.genfromtxt('1outP0_R_woInt_10m6_wInt_mN_y10_rc2_stable.dat'))
data.append(np.genfromtxt('1outP0_R_woInt_10m6_wInt_mN_y10_rc4_stable.dat'))
data.append(np.genfromtxt('1outP0_R_woInt_10m6_wInt_mN_y10_rc6_stable.dat'))
data.append(np.genfromtxt('1outP0_R_woInt_10m6_wInt_mN_y10_rc8_stable.dat'))
#data.append(np.genfromtxt('1outP0_R.dat'))#_woInt_10m6_wInt_mN_y10_rc4.dat'))


fig, ax = plt.subplots()


ax.plot(data[0][:,0], data[0][:,1], label = "$r_c=1$" , color= 'blue')#'#de0c62')	
ax.plot(data[1][:,0], data[1][:,1], label = "$r_c=2$" , color= '#438900')#'#5fa052')	
ax.plot(data[2][:,0], data[2][:,1], label = "$r_c=4$" , color= '#00a9df')	
ax.plot(data[3][:,0], data[3][:,1], label = "$r_c=6$" , color= 'violet')#'#ff8d0d')	
ax.plot(data[4][:,0], data[4][:,1], label = "$r_c=8$" , color= 'purple')#'#071cf4')


ax.text(25,0.008,
         r'$y=10$' , fontsize=15 ,
         bbox={'facecolor': 'white', 'alpha': 0.3, 'pad': 6})


ax.set_xlabel(r'$R^{\prime}$', fontsize=18)
ax.set_xlim(1,60)#1E+2)
ax.set_xscale('log')

ax.set_ylabel(r'$P^{\prime}_{c}$ ', fontsize=18, usetex=True)
ax.set_ylim(1E-11,0.1)#-0.39,0.29)#-0.0999,-0.0996)#0)#1.3)   09965
ax.set_yscale('log')

ax.legend(loc='lower left', fontsize=16, frameon=True)#, mode = "expand", ncol = 3)

ax.minorticks_on()



plt.savefig("P0_R_fermi_woInt_10m6_wInt_mN_y10_stable1.pdf",dpi=300,bbox_inches='tight',pad_inches=0.015,transparent=False)