
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 10:35:50 2021

@author: Benutzer01
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

figw=5#14#5.
figh=4#8#5.

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
#data.append(np.genfromtxt('1outP_Pc1_rho05_rc01.dat'))
#data.append(np.genfromtxt('1outP_Pc1_rho05_rc02.dat'))
#data.append(np.genfromtxt('1outP_Pc1_rho05_rc025.dat'))
#data.append(np.genfromtxt('1outP_Pc1_rho05_rc03.dat'))
data.append(np.genfromtxt('1outM_R_rc005_rho05_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc01_rho05_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc015_rho05_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc02_rho05_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc025_rho05_pc_var_sound.dat'))
#data.append(np.genfromtxt('1outM_R_rc01_rho05_pc_var_sound.dat'))
#data.append(np.genfromtxt('1outM_R.dat'))

fig, ax = plt.subplots( figsize=(figw,figh),dpi=400)

xdat = np.linspace(0,2,1000)
ax.plot(xdat, 3/8*xdat, "--", color= 'black',linewidth=1.5)
ax.plot(xdat, 4/9*xdat, "-", color= 'black')



ax.plot(data[0][:,0], data[0][:,1], label = "$r_c=0.05$" , color= 'red')#'#de0c62')	
ax.plot(data[1][:,0], data[1][:,1], label = "$r_c=0.1$" , color= '#5fa052')	
ax.plot(data[2][:,0], data[2][:,1], label = "$r_c=0.15$" , color= '#00a9df')	
ax.plot(data[3][:,0], data[3][:,1], label = "$r_c=0.2 $ " , color= 'blue')#'#ff8d0d')	
ax.plot(data[4][:,0], data[4][:,1], label = "$r_c=0.25 $" , color= 'orange')#'#071cf4')
#ax.plot(data[5][:,0], data[5][:,1], label = "$r_c=1$ sound" , color= '#d76eee')
#ax.plot(data[6][:,0], data[5][:,1], label = "$r_c=1$ sound" , color= '#d76eee')


#ax.text(0.83,2E+1, #0.22 0.05, 0.12, #'$\epsilon  = 0$' "\n" 
#         r'$\epsilon_{core} = 10^{-4} $' "\n" 
#         r'$m_{shell} = 10^{-6}m_N$' , fontsize=15  ) #"\n" )
#         r'$a_1 = 0$' , fontsize=15 )
#         bbox={'facecolor': 'white', 'alpha': 0.3, 'pad': 6})

#ax.text(30.83,0.5, #0.22 0.05, 0.12, #'$\epsilon  = 0$' "\n" 
#         r'$\epsilon_{core} = 10^{-4} $' "\n" 
#         r'$m_{shell} = 10^{-6}m_N$' , fontsize=15  ) #"\n" )
#         r'$y=10$' , fontsize=15 ,
#         bbox={'facecolor': 'white', 'alpha': 0.3, 'pad': 6})
ax.text(0.355,0.015, 
         r'$\rho_1 / \rho_0 = 0.5$' , fontsize=14,
         bbox={'facecolor': 'white', 'alpha': 0.3, 'pad': 6})

#ax.text(0.03,0.077, 
#         r'$\rho_1 / \rho_0 = 0.5$' , fontsize=14 )

ax.set_xlabel(r'$R^{\prime}$', fontsize=18)
ax.set_xlim(0,0.5) #1E-1,1E+2)
#ax.set_xscale('log')

ax.set_ylabel(r'$M^{\prime}$ ', fontsize=18, usetex=True)
ax.set_ylim(0,0.185)#1E-11,1E+1)#-0.39,0.29)#-0.0999,-0.0996)#0)#1.3)   09965
#ax.set_yscale('log')

#ax.set_title(r"The effective potential of Schwarzschild spacetime")

ax.legend(loc='upper left', fontsize=14, frameon=True)#, mode = "expand", ncol = 3)

ax.minorticks_on()



plt.savefig("M_rc_diff_rho_05_Pc_var_causality_sound.pdf",dpi=300,bbox_inches='tight',pad_inches=0.015,transparent=False)