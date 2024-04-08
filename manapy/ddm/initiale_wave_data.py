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




#TODO initialisation
@njit(fastmath=True)
def initialisation_wave(w_c, center, choix_wave, ho, grav, k1, k2, B01, B02, tol, wave):

    nbelements = len(center)
    u0=0
    v0=0    
    B0 = np.sqrt(B01**2 + B02**2)
    k  = np.sqrt(k1**2  + k2**2)
    c0 = np.sqrt(ho*grav)
    # Produit scalaire de (k1,k2) et (B01, B02)
    SS = k1*B01 + k2*B02
    theta = np.arccos(SS/(B0*k))
    w1 = k*np.sqrt(c0*c0 + B0*np.cos(theta)) 
    w2 = k*B0*np.fabs(np.cos(theta)) 

    
    if wave == "MG" :
        h_hot  = ho*k*k
        u_hot  = -k1*w1
        v_hot  = -k2*w1
        B1_hot = -k1*w2
        B2_hot = -k2*w2
        
    if wave == "AL" :
        h_hot  = 0.0
        u_hot  = k2
        v_hot  = -k1
        B1_hot = k2
        B2_hot = -k1
                
    if choix_wave == 1 : 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            if 0.5 < xcent and xcent < 1 :
            	  tol1 = tol             
            else  :
            	  tol1 = 0.0                         

            
            hint  = ho  + tol1*h_hot *np.cos(k1*xcent + k2*ycent)
            uint  = u0  + tol1*u_hot * np.cos(k1*xcent + k2*ycent)
            vint  = v0  + tol1*v_hot * np.cos(k1*xcent + k2*ycent)
            B1int = B01 + tol1*B1_hot* np.cos(k1*xcent + k2*ycent)
            B2int = B02 + tol1*B2_hot* np.cos(k1*xcent + k2*ycent)
      
            w_c.h[i]  = hint
            w_c.hu[i] =  hint*uint 
            w_c.hv[i] =  hint*vint 
            w_c.hB1[i] = hint*B1int 
            w_c.hB2[i] = hint*B2int            

    elif choix_wave == 2 : 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]

            if  xcent >=49 and xcent <= 51 :
            	  tol1 = tol             
            else  :
            	  tol1 = 0.0   
            hint  = tol1*h_hot *np.cos(k1*xcent + k2*ycent)
            uint  = tol1*u_hot * np.cos(k1*xcent + k2*ycent)
            vint  = tol1*v_hot * np.cos(k1*xcent + k2*ycent)
            B1int = tol1*B1_hot* np.cos(k1*xcent + k2*ycent)
            B2int = tol1*B2_hot* np.cos(k1*xcent + k2*ycent)
      
            AA = 0.0*np.sin(2*np.pi*(xcent-ycent))


            if xcent >=49 and xcent <= 51 : 
                     aa =100 + 0.1
            else  :                
                     aa = 100.

            uu = 2.
            bb = 5.#*0.0
            w_c.h[i]  = aa + hint
            w_c.hu[i] = uu + AA + hint*uint
            w_c.hv[i] = uu + AA + hint*vint 
            w_c.hB1[i] = bb + 2*AA + hint*B1int 
            w_c.hB2[i] = bb + 2*AA + hint*B2int  

            w_c.hP[i] = 0.   
            w_c.Ch[i] = 0. 


    return w_c


