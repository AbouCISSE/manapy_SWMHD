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
def initialisation_source_terme(w_c, center, choix, ho, grav, k1, k2, B01, B02, tol, wave):

    nbelements = len(center)




    if choix == 1: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            w_c.h[i]  = 2. + np.sin(2*np.pi*xcent) 
            w_c.hu[i] = 2. + np.sin(2*np.pi*xcent) 
            w_c.hv[i] = 2. + np.sin(2*np.pi*xcent) 
            w_c.hB1[i] = 1.
            w_c.hB2[i] = 4 + np.sin(2*np.pi*xcent) 
                            
   

    elif choix == 2: #(c-propriété) [-10 10]X[-10 10]
       
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            a1 = 0.75*np.exp(-0.5*((xcent+4)**2 +(ycent-5)**2))+ 0.7*np.exp(-0.5*((xcent+2.5)**2 +(ycent-2.5)**2))
            a2 = 0.65*np.exp(-0.3*(xcent**2 +ycent**2))+ 0.6*np.exp(-0.4*((xcent-3)**2 +(ycent+2)**2))  
            a3 = 0.55*np.exp(-0.7*((xcent-5)**2 +(ycent+4)**2)) 
            
            w_c.Z[i] = a1 + a2 + a3
            w_c.h[i] = ho - w_c.Z[i]
            w_c.hu[i] = 0.
            w_c.hv[i] = 0.
            w_c.hB1[i] = 0.
            w_c.hB2[i] = 0.  
            w_c.hP[i] = 0.


    elif choix == 3: # (c-propriété) h gaussienne  [-0,2]X[0,1]

        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            if xcent > 0.05 and xcent < 0.15:
                
                w_c.h[i] = ho + 0.01
            
            else:
                  w_c.h[i] = ho
                  
            w_c.Z[i] = 0.8 * np.exp(-5*(xcent - 0.9)**2 - 50* (ycent - 0.5)**2) 
            w_c.h[i] =  w_c.h[i] - w_c.Z[i]
            w_c.hu[i] = 0.
            w_c.hv[i] = 0.
            w_c.hB1[i] = 0.
            w_c.hB2[i] = 0.  
            w_c.hP[i] = 0.
            
    elif choix == 4:
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]

            w_c.Z[i] = (10*np.exp(-1*((xcent-2)**2+(ycent-6)**2)) + 15*np.exp(-1*((xcent-2.5)**2+(ycent-2.5)**2)))/20
            w_c.Z[i] += (12*np.exp(-1*((xcent-5)**2+(ycent-5)**2)) + 6*np.exp(-2*((xcent-7.5)**2+(ycent-7.5)**2)))/20
            w_c.Z[i] += (16*np.exp(-1*((xcent-7.5)**2+(ycent-2)**2)))/20

            w_c.h[i] = 1 - w_c.Z[i]
            w_c.hu[i] = 0.
            w_c.hv[i] = 0.
            w_c.hB1[i] = 1
            w_c.hB2[i] = 1.  
            w_c.hP[i] = 0.

    elif choix == 5: #h gaussienne (c-propriété)

        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]

            w_c.Z[i] = 0.8 * np.exp(-5*(xcent - 1)**2 - 50* (ycent - 0.5)**2) 

            
            if 0.05 < xcent and xcent < 0.15 :
            	   w_c.h[i] = 1. - w_c.Z[i] # 1.01
            else :
            	   w_c.h[i] = 1 - w_c.Z[i]

            w_c.hu[i] = 0.
            w_c.hv[i] = 0.
            w_c.hB1[i] = 0.
            w_c.hB2[i] = 0.
            w_c.hP[i] = 0.
            
    elif choix == 6: 

        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            g = 1
            s1 = np.sin(2*np.pi*xcent)
            s2 = s1 + 2
            s3 = 2*g

            w_c.Z[i] = 0# s1 - 1/(s3*pow(s2,2)) + 2

            w_c.h[i] = 2 + s1
            w_c.hu[i] = 2 + s1
            w_c.hv[i] = 2 + s1
            w_c.hB1[i] = 1
            w_c.hB2[i] = 4 + 2*s1
            w_c.hP[i] = 0.

    elif choix == 7: #h gaussienne (c-propriété)

        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]

            w_c.Z[i] = 0.8 * np.exp(-5*(xcent - 1)**2 - 50* (ycent - 0.5)**2) 

            
            if 0.05 < xcent and xcent < 0.5 :
            	   hint = 1.01 - w_c.Z[i] # 1.01
            else :
            	   hint = 1 - w_c.Z[i]

            w_c.h[i] = hint
            w_c.hu[i] = 0.
            w_c.hv[i] = 0.
            w_c.hB1[i] = 0.*1.002*hint
            w_c.hB2[i] = 0.*0.1*hint
            w_c.hP[i] = 0.

    return w_c


