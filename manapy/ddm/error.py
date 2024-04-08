#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 19:04:50 2020

@author: kissami
"""
import meshio

import numpy as np
from numba import njit, jit
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()




@njit(fastmath=True)
def global_cpsi(w_c, cfl, normal, mesure, volume, faceid, grav, globav):
    nbelement =  len(faceid)
    dt_c = 1e6
    kl=10
    c1 = []
    c2 = 0

    for i in range(nbelement):
        velson = np.sqrt(grav*w_c[i].h)
        
 
        if globav == "on" :
            for j in range(3):
                norm = normal[faceid[i][j]]
                mesn = mesure[faceid[i][j]]
                u_n = w_c.hu[i]/w_c.h[i]*norm[0] + w_c.hv[i]/w_c.h[i]*norm[1]
                u_n = u_n/mesn
                B_n = w_c.hB1[i]/w_c.h[i]*norm[0] + w_c.hB2[i]/w_c.h[i]*norm[1]
                B_n = B_n/mesn
                wb  = np.sqrt(B_n**2 + grav*w_c[i].h)
                lam1 = np.fabs(u_n - wb)





                    
def norml2(w_c, vol, cells, nodes, w_x, w_y, ET):


    cells  = np.array(cells)
    nbelement = len(w_c)
    divhB = np.zeros(len(w_c))
    ETL   = np.zeros(len(w_c))
    eng=0.0


    for i in range(nbelement):

        divhB[i] = (w_x["hB1"][i] + w_y["hB2"][i])*vol[i]
        ETL[i]   = ET[i]*vol[i]
        eng += (pow(ET[i],2))*vol[i]
        

        if globav == "on" :
            for j in range(3):
                norm = normal[faceid[i][j]]
                mesn = mesure[faceid[i][j]]
   
                B_n = (w_c.hB1[i]*norm[0] + w_c.hB2[i]*norm[1])/mesn
                lam1 = np.fabs(u_n - wb)




    dErrorL2 = np.linalg.norm(divhB,ord=2)
    div = dErrorL2
    #eng = np.linalg.norm(ETL,ord=2)
    
    return div, np.sqrt(eng)


def Total_Energy(w_c, grav):

    nbelement = len(w_c)
    Et = np.zeros(nbelement) # Total_Energy
    Ec = np.zeros(nbelement) # Kinetic Energy
    Ep = np.zeros(nbelement) # Potential Energy
    Em = np.zeros(nbelement) # Magnetic energy
    
    for i in range(nbelement):

        hc  = w_c["h"][i]
        uc  = w_c["hu"][i]/w_c["h"][i]
        vc  = w_c["hv"][i]/w_c["h"][i]
        B1c = w_c["hB1"][i]/w_c["h"][i]
        B2c = w_c["hB2"][i]/w_c["h"][i]
        bc  = w_c["Z"][i]
        Ec[i] = 0.5*hc*(pow(uc , 2)  + pow(vc , 2))
        Em[i] = 0.5*hc*(pow(B1c , 2) + pow(B2c , 2))
        Ep[i] = 0.5*grav*pow(hc , 2) + grav*bc*hc
        Et[i] = Ec[i] + Em[i] + Ep[i]
  
    return Et

