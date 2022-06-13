
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
data.append(np.genfromtxt('1outM_R_rc02_rho01_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc02_rho02_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc02_rho03_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc02_rho05_pc_var_sound.dat'))
data.append(np.genfromtxt('1outM_R_rc02_rho08_pc_var_sound.dat'))
#data.append(np.genfromtxt('1outM_R.dat'))

fig, ax = plt.subplots( figsize=(figw,figh),dpi=400)

xdat = np.linspace(0,2,1000)
ax.plot(xdat, 3/8*xdat, "--", color= 'black',linewidth=1.5)
ax.plot(xdat, 4/9*xdat, "-", color= 'black')


ax.plot(data[0][:,0], data[0][:,1], label = r"$ {\rho}_{1} / {\rho}_0 =0.1$" , color= 'red')#'#de0c62')	
ax.plot(data[1][:,0], data[1][:,1], label = r"${\rho}_{1} / {\rho}_0 =0.2$" , color= 'orange')#'#5fa052')	
ax.plot(data[2][:,0], data[2][:,1], label = r"${\rho}_{1} / {\rho}_0 =0.3$" , color= '#00a9df')	
ax.plot(data[3][:,0], data[3][:,1], label = r"${\rho}_{1} / {\rho}_0 =0.5$" , color= 'purple')#'#ff8d0d')	
ax.plot(data[4][:,0], data[4][:,1], label = r"${\rho}_{1} / {\rho}_0 =0.8$" , color= '#029f06')#'#071cf4')
#ax.plot(data[5][:,0], data[5][:,1], label = "$r_c=0.28$" , color= '#d76eee')

#ax.plot(data[6][:,0], data[6][:,1], label = "$r_c=7$" , color= '#6eeed1')
#ax.plot(data[7][:,0], data[7][:,1], label = "$r_c=8$" , color= '#b80808')
#ax.plot(data[8][:,0], data[8][:,1], label = "$r_c=9$" , color= '#ff8b00')
#ax.plot(data[9][:,0], data[9][:,1], label = "$r_c=10$" , color= '#8b04f5')
#ax.plot(data[10][:,0], data[10][:,1], label = "$y100$" , color= '#ff8b00')

#ax.vlines(x=0.3475,  ymin=-5, ymax=-0.08613, ls=':', color='#5fa052')  # cusps
#ax.vlines(x=0.60,  ymin=-5, ymax=-0.08613, ls=':', color='#5fa052')
#ax.vlines(x=0.7925,  ymin=-5, ymax=-0.08613, ls=':', color='#5fa052')
#ax.hlines(xmin=0.3475, xmax=0.7925, y=-0.08613, ls=':', color='#5fa052')

#ax.vlines(x=0.4475,  ymin=-5, ymax=-0.08898, ls=':', color=clr[2])#'#de0c62')  # maxima
#ax.vlines(x=0.7325,  ymin=-5, ymax=-0.08831, ls=':', color=clr[2])#'#de0c62')

#ax.vlines(x=0.315,  ymin=-5, ymax=-0.06692, ls=':', color=clr[2])    #in out
#ax.vlines(x=0.508,  ymin=-5, ymax=-0.06692, ls=':', color=clr[2])
#ax.hlines(xmin=0.315, xmax=0.508, y=-0.06692, ls=':', color=clr[2])

#ax.vlines(x=0.740,  ymin=-5, ymax=-0.06489, ls=':', color=clr[2])   # in out
#ax.vlines(x=0.890,  ymin=-5, ymax=-0.06489, ls=':', color=clr[2])
#ax.hlines(xmin=0.740, xmax=0.890, y=-0.06489, ls=':', color=clr[2])
    
#ax.text(0.338, -0.0859, r'$r_{\mathrm{cusp}}$', fontsize=14)
#ax.text(0.58, -0.0859, r'$r_{\mathrm{cusp}}$', fontsize=14)
#ax.text(0.37, -0.0893, r'$r_{\mathrm{max}}$', fontsize=14)
#ax.text(0.65, -0.0887, r'$r_{\mathrm{max}}$', fontsize=14)
#plt.rcParams['mathtext.rm'] = 'Modern'#'stix:italic' 

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
ax.text(0.79,0.04, 
         r'$r_c = 0.2$' , fontsize=14,
         bbox={'facecolor': 'white', 'alpha': 0.3, 'pad': 6})

#ax.text(0.05,0.195, 
 #        r'$r_c = 0.2$' , fontsize=14 )

ax.set_xlabel(r'$R^{\prime}$', fontsize=18)
ax.set_xlim(0,1) #1E-1,1E+2)
#ax.set_xscale('log')

ax.set_ylabel(r'$M^{\prime}$ ', fontsize=18, usetex=True)
ax.set_ylim(0,0.5)#1E-11,1E+1)#-0.39,0.29)#-0.0999,-0.0996)#0)#1.3)   09965
#ax.set_yscale('log')

#ax.set_title(r"The effective potential of Schwarzschild spacetime")

ax.legend(loc='upper left', fontsize=14, frameon=True)#, mode = "expand", ncol = 3)

ax.minorticks_on()



plt.savefig("M_rc02_rho_diff_Pc_var_causality_sound.pdf",dpi=300,bbox_inches='tight',pad_inches=0.015,transparent=False)