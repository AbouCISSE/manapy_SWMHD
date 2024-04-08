#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 20:51:37 2020

@author: kissami
"""
import timeit
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
plt.rcParams['text.usetex'] = True


E1 = np.loadtxt("k0.txt", dtype='float', delimiter=',').T
E2 = np.loadtxt("k2.txt", dtype='float', delimiter=',').T
E3 = np.loadtxt("k20.txt", dtype='float', delimiter=',').T
E4 = np.loadtxt("k60.txt", dtype='float', delimiter=',').T


E1B1 = E1[0] ; E1B2 = E1[1] ; E1h  = E1[2] ; E1u  = E1[4] ; E1v  = E1[5] ; E1x  = E1[9] ; E1y  = E1[10]  
E2B1 = E2[0] ; E2B2 = E2[1] ; E2h  = E2[2] ; E2u  = E2[4] ; E2v  = E2[5] ; E2x  = E2[9] ; E2y  = E2[10]  
E3B1 = E3[0] ; E3B2 = E3[1] ; E3h  = E3[2] ; E3u  = E3[4] ; E3v  = E3[5] ; E3x  = E3[9] ; E3y  = E3[10]  
E4B1 = E4[0] ; E4B2 = E4[1] ; E4h  = E4[2] ; E4u  = E4[4] ; E4v  = E4[5] ; E4x  = E4[9] ; E4y  = E4[10]  


fig, ax = plt.subplots(1, 1,)# figsize=(7,7))
#ax.plot(E1x, E1h, "k",linewidth=2, markersize=6, label='t=0s')
ax.plot(E2x, E2h, "k*:",linewidth=2, markersize=1, label='t=2s')
ax.plot(E3x, E3h, "r",linewidth=2, markersize=1, label='t=20s')
ax.plot(E4x, E4h, "b",linewidth=2, markersize=6, label='t=60s')


#ax.yaxis.set_major_formatter(FormatStrFormatter('%0.15g'))

#plt.yscale('log')
plt.title(r'case $\mathbf{k} \perp (Ox)$', fontsize=20)
plt.xlabel(r"$x$", fontsize=20)
plt.ylabel(r"$h$", fontsize=20)
#plt.ylim([96,104])
plt.xlim([0, 10])
#plt.ylabel(r'$E_T(t)$')
ax.legend(fontsize=13, loc="best")
left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
#ax2 = fig.add_axes([left, bottom, width, height])
#ax2.plot(E3x, E3h, "g",linewidth=2, markersize=6)
#ax2.plot(E4x, E4h, "b-.",linewidth=2, markersize=6)
#ax2.tick_params(top=False, bottom=False, left=False, right=False,
#                labelleft=False, labelbottom=False)
#ax2.yaxis.set_major_formatter(FormatStrFormatter('%0.15g'))
plt.savefig("hky1.png") 
plt.show()
