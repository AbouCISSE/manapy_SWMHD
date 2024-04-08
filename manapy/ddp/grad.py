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
def cell_gradient2(w_c, w_ghost, w_halo, w_haloghost, centerc, vertexn, nodeid, cellnid, halonid, halocenterg, namen, centerg, centerh, w_x, w_y,shift, periodicn):

    nbelement = len(centerc)
    center = np.zeros(3)
    for i in range(nbelement):
        i_xx = 0.
        i_yy = 0.
        i_xy = 0.
        j_xh = 0.
        j_yh = 0.

        j_xhu = 0.
        j_yhu = 0.
        j_xhv = 0.
        j_yhv = 0.
        j_xhB1 = 0.
        j_yhB1 = 0.
        j_xhB2 = 0.
        j_yhB2 = 0.
        j_xz = 0.
        j_yz = 0.

        for j in range(len(cellnid[i])):
            cell = cellnid[i][j]
            if cell != -1:
                j_x = centerc[cell][0] - centerc[i][0]
                j_y = centerc[cell][1] - centerc[i][1]
                i_xx += pow(j_x, 2)
                i_yy += pow(j_y, 2)
                i_xy += (j_x * j_y)

                j_xh += (j_x * (w_c[cell].h - w_c[i].h))
                j_yh += (j_y * (w_c[cell].h - w_c[i].h))

                j_xhu += (j_x * (w_c[cell].hu - w_c[i].hu))
                j_yhu += (j_y * (w_c[cell].hu - w_c[i].hu))

                j_xhv += (j_x * (w_c[cell].hv - w_c[i].hv))
                j_yhv += (j_y * (w_c[cell].hv - w_c[i].hv))

                j_xhB1 += (j_x * (w_c[cell].hB1 - w_c[i].hB1))
                j_yhB1 += (j_y * (w_c[cell].hB1 - w_c[i].hB1))
                
                j_xhB2 += (j_x * (w_c[cell].hB2 - w_c[i].hB2))
                j_yhB2 += (j_y * (w_c[cell].hB2 - w_c[i].hB2))
                
                j_xz += (j_x * (w_c[cell].Z - w_c[i].Z))
                j_yz += (j_y * (w_c[cell].Z - w_c[i].Z))



        for j in range(len(periodicn[i])):
            cell = periodicn[i][j]
            
            if cell != -1:
                #print("I'm okay")
                center[:] = centerc[cell][:]                
                j_x = center[0]  - centerc[i][0] + shift[cell][0]
                j_y = center[1]  - centerc[i][1] + shift[cell][1]
                i_xx += j_x*j_x
                i_yy += j_y*j_y
                i_xy += (j_x * j_y)
                
                
                j_xh += (j_x * (w_c[cell].h - w_c[i].h))
                j_yh += (j_y * (w_c[cell].h - w_c[i].h))

                j_xhu += (j_x * (w_c[cell].hu - w_c[i].hu))
                j_yhu += (j_y * (w_c[cell].hu - w_c[i].hu))

                j_xhv += (j_x * (w_c[cell].hv - w_c[i].hv))
                j_yhv += (j_y * (w_c[cell].hv - w_c[i].hv))

                j_xhB1 += (j_x * (w_c[cell].hB1 - w_c[i].hB1))
                j_yhB1 += (j_y * (w_c[cell].hB1 - w_c[i].hB1))
                
                j_xhB2 += (j_x * (w_c[cell].hB2 - w_c[i].hB2))
                j_yhB2 += (j_y * (w_c[cell].hB2 - w_c[i].hB2))
                
                j_xz += (j_x * (w_c[cell].Z - w_c[i].Z))
                j_yz += (j_y * (w_c[cell].Z - w_c[i].Z))
                
  
        
        for k in range(3):
            nod = nodeid[i][k]
            if vertexn[nod][3] <= 4:
                for j in range(len(centerg[nod])):
                    cell = int(centerg[nod][j][-1])
                    if cell != -1:
                        center[:] = centerg[nod][j][0:3]
                        j_x = center[0] - centerc[i][0]
                        j_y = center[1] - centerc[i][1]
                        
                        i_xx += j_x*j_x
                        i_yy += j_y*j_y
                        i_xy += (j_x * j_y)
                        j_xh += (j_x * (w_ghost[cell].h - w_c[i].h))
                        j_yh += (j_y * (w_ghost[cell].h - w_c[i].h))

                        j_xhu += (j_x * (w_ghost[cell].hu - w_c[i].hu))
                        j_yhu += (j_y * (w_ghost[cell].hu - w_c[i].hu))

                        j_xhv += (j_x * (w_ghost[cell].hv - w_c[i].hv))
                        j_yhv += (j_y * (w_ghost[cell].hv - w_c[i].hv))

                        j_xhB1 += (j_x * (w_ghost[cell].hB1 - w_c[i].hB1))
                        j_yhB1 += (j_y * (w_ghost[cell].hB1 - w_c[i].hB1))
                     
                        j_xhB2 += (j_x * (w_ghost[cell].hB2 - w_c[i].hB2))
                        j_yhB2 += (j_y * (w_ghost[cell].hB2 - w_c[i].hB2))

                        j_xz += (j_x * (w_ghost[cell].Z - w_c[i].Z))
                        j_yz += (j_y * (w_ghost[cell].Z - w_c[i].Z))                 
                    #    j_xw += (j_x * (w_ghost[cell] - w_c[i] ))
                     #   j_yw += (j_y * (w_ghost[cell] - w_c[i] ))

                        
            if namen[nod] == 10 :
                for j in range(len(halocenterg[nod])):
                    #-3 the index of global face
                    cell = int(halocenterg[nod][j][-1])
                    if cell != -1:
                        center[:] = halocenterg[nod][j][:]
                        j_x = center[0] - centerc[i][0]
                        j_y = center[1] - centerc[i][1]
                        
                        i_xx += j_x*j_x
                        i_yy += j_y*j_y
                        i_xy += (j_x * j_y)
        
                        j_xh += (j_x * (w_haloghost[cell].h - w_c[i].h ))
                        j_yh += (j_y * (w_haloghost[cell].h - w_c[i].h ))

                        j_xhu += (j_x * (w_haloghost[cell].hu - w_c[i].hu ))
                        j_yhu += (j_y * (w_haloghost[cell].hu - w_c[i].hu ))
  
                        j_xhv += (j_x * (w_haloghost[cell].hv - w_c[i].hv ))
                        j_yhv += (j_y * (w_haloghost[cell].hv - w_c[i].hv ))

                        j_xhB1 += (j_x * (w_haloghost[cell].hB1 - w_c[i].hB1 ))
                        j_yhB1 += (j_y * (w_haloghost[cell].hB1 - w_c[i].hB1 ))

                        j_xhB2 += (j_x * (w_haloghost[cell].hB2 - w_c[i].hB2 ))
                        j_yhB2 += (j_y * (w_haloghost[cell].hB2 - w_c[i].hB2 ))

                        j_xz += (j_x * (w_haloghost[cell].Z - w_c[i].Z ))
                        j_yz += (j_y * (w_haloghost[cell].Z - w_c[i].Z ))

        dia = i_xx * i_yy - pow(i_xy, 2)

        w_x[i].h = (i_yy * j_xh - i_xy * j_yh) / dia
        w_y[i].h = (i_xx * j_yh - i_xy * j_xh) / dia

        w_x[i].hu = (i_yy * j_xhu - i_xy * j_yhu) / dia
        w_y[i].hu = (i_xx * j_yhu - i_xy * j_xhu) / dia

        w_x[i].hv = (i_yy * j_xhv - i_xy * j_yhv) / dia
        w_y[i].hv = (i_xx * j_yhv - i_xy * j_xhv) / dia

        w_x[i].hB1 = (i_yy * j_xhB1 - i_xy * j_yhB1) / dia
        w_y[i].hB1 = (i_xx * j_yhB1 - i_xy * j_xhB1) / dia

        w_x[i].hB2 = (i_yy * j_xhB2 - i_xy * j_yhB2) / dia
        w_y[i].hB2 = (i_xx * j_yhB2 - i_xy * j_xhB2) / dia

        w_x[i].Z = (i_yy * j_xz - i_xy * j_yz) / dia
        w_y[i].Z = (i_xx * j_yz - i_xy * j_xz) / dia

    return w_x, w_y














def centertovertex(w_c, w_ghost, w_halo, w_haloghost, centerc, centerh, cellidn, periodicn, haloidn, vertexn,  
                   namen, centergn, halocentergn, R_x, R_y, lambda_x,lambda_y, number, shift, nbproc,  w_n):
   
    from numpy import zeros
    nbnode = len(vertexn)
    center = zeros(3)
    
    for i in range(nbnode):
        for j in range(cellidn[i][-1]):
            cell = cellidn[i][j]
            center[:] = centerc[cell][:]
           
            xdiff = center[0] - vertexn[i][0]
            ydiff = center[1] - vertexn[i][1]
            alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
       
            w_n[i]  += alpha * w_c[cell]
                
        # if vertexn[i][3] == 1 or vertexn[i][3] == 2 or vertexn[i][3] == 3 or vertexn[i][3] == 4:
        if centergn[i][0][2] != -1: 
            for j in range(len(centergn[i])):
                cell = int(centergn[i][j][-1])
                if cell != -1:
                    #center[:] = centergn[i][j][:]
                    
                    #xdiff = center[0] - vertexn[i][0]
                    #ydiff = center[1] - vertexn[i][1]
                    xdiff = centergn[i][j][0] - vertexn[i][0]
                    ydiff = centergn[i][j][1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                    
                    w_n[i]  += alpha * w_ghost[cell]
                    
        if vertexn[i][3] == 5 or vertexn[i][3] == 6 :
            for j in range(len(periodicn[i])):
                cell = periodicn[i][j]
                if cell != -1:
                    center[:] = centerc[cell][:] 
                    
                    xdiff = center[0] + shift[cell][0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                    
                    w_n[i]  += alpha * w_c[cell]
                    
        elif vertexn[i][3] == 7 or vertexn[i][3] == 8:
            for j in range(len(periodicn[i])):
                cell = periodicn[i][j]
                if cell != -1:
                    center[:] = centerc[cell][:] 
                    
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] + shift[cell][1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                    
                    w_n[i]  += alpha * w_c[cell]
                    
        if namen[i] == 10:
           
            for j in range(len(halocentergn[i])):
                cell = int(halocentergn[i][j][-1])
                if cell != -1:
                    center[:] = halocentergn[i][j][:]
                  
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                 
                    w_n[i]  += alpha * w_haloghost[cell]
            
            if haloidn[i][-1] > 0 :
                for j in range(haloidn[i][-1]):
                    cell = haloidn[i][j]
                    center[:] = centerh[cell][:]
                  
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                 
                    w_n[i]  += alpha * w_halo[cell]































def cell_gradient(w_c, w_ghost, w_halo, w_haloghost, centerc, cellnid, halonid,
                  nodecid, periodicn, namen, centerg, halocenterg, vertexn, centerh, shift,
                  nbproc, w_x, w_y):
                          
    
    from numpy import zeros
    center = zeros(3)
    nbelement = len(w_c)
    
    for i in range(nbelement):
        i_xx  = 0.;  i_yy  = 0.; i_xy = 0.
        j_xw = 0.;  j_yw = 0.

        for j in range(cellnid[i][-1]):
            cell = cellnid[i][j]
            j_x = centerc[cell][0] - centerc[i][0]
            j_y = centerc[cell][1] - centerc[i][1]
            i_xx += j_x*j_x
            i_yy += j_y*j_y
            i_xy += (j_x * j_y)

            j_xw += (j_x * (w_c[cell] - w_c[i] ))
            j_yw += (j_y * (w_c[cell] - w_c[i] ))
            
        for j in range(len(periodicn[i])):
            cell = periodicn[i][j]
            if cell != -1:
                center[:] = centerc[cell][:]
                j_x = center[0] + shift[cell][0] - centerc[i][0]
                j_y = center[1] + shift[cell][1] - centerc[i][1]
                
                i_xx += j_x*j_x
                i_yy += j_y*j_y
                i_xy += (j_x * j_y)
                
                j_xw += (j_x * (w_c[cell] - w_c[i] ))
                j_yw += (j_y * (w_c[cell] - w_c[i] ))
                
        # if nbproc > 1:      
        for j in range(halonid[i][-1]):
            cell = halonid[i][j]
            j_x = centerh[cell][0] - centerc[i][0]
            j_y = centerh[cell][1] - centerc[i][1]
            
            i_xx += j_x*j_x
            i_yy += j_y*j_y
            i_xy += (j_x * j_y)
            
            j_xw += (j_x * (w_halo[cell]  - w_c[i] ))
            j_yw += (j_y * (w_halo[cell]  - w_c[i] ))
                
                
        # #TODO verify ghost center
        for k in range(3):
            nod = nodecid[i][k]
            if vertexn[nod][3] <= 4:
                for j in range(len(centerg[nod])):
                    cell = int(centerg[nod][j][-1])
                    if cell != -1:
                        center[:] = centerg[nod][j][0:3]
                        j_x = center[0] - centerc[i][0]
                        j_y = center[1] - centerc[i][1]
                        
                        i_xx += j_x*j_x
                        i_yy += j_y*j_y
                        i_xy += (j_x * j_y)
                        
                        j_xw += (j_x * (w_ghost[cell] - w_c[i] ))
                        j_yw += (j_y * (w_ghost[cell] - w_c[i] ))
                        
            if namen[nod] == 10 :
                for j in range(len(halocenterg[nod])):
                    #-3 the index of global face
                    cell = int(halocenterg[nod][j][-1])
                    if cell != -1:
                        center[:] = halocenterg[nod][j][:]
                        j_x = center[0] - centerc[i][0]
                        j_y = center[1] - centerc[i][1]
                        
                        i_xx += j_x*j_x
                        i_yy += j_y*j_y
                        i_xy += (j_x * j_y)
        
                        j_xw += (j_x * (w_haloghost[cell] - w_c[i] ))
                        j_yw += (j_y * (w_haloghost[cell] - w_c[i] ))

        dia = i_xx * i_yy - i_xy*i_xy

        w_x[i]  = (i_yy * j_xw - i_xy * j_yw) / dia
        w_y[i]  = (i_xx * j_yw - i_xy * j_xw) / dia
        
    return w_x, w_y

@njit(fastmath=True)
def variablesNP(centerc, cellidn, haloidn, vertexn, namen, centerg, centerh):
    
    nbnode = len(vertexn)

    I_xx = np.zeros(nbnode)
    I_yy = np.zeros(nbnode)
    I_xy = np.zeros(nbnode)
    R_x = np.zeros(nbnode)
    R_y = np.zeros(nbnode)
    lambda_x = np.zeros(nbnode)
    lambda_y = np.zeros(nbnode)
    number = np.zeros(nbnode)


    for i in range(nbnode):
        for j in range(cellidn[i][-1]):
            center = centerc[cellidn[i][j]]
            Rx = center[0] - vertexn[i][0]
            Ry = center[1] - vertexn[i][1]
            I_xx[i] += (Rx * Rx)
            I_yy[i] += (Ry * Ry)
            I_xy[i] += (Rx * Ry)
            R_x[i] += Rx
            R_y[i] += Ry
            number[i] = number[i] + 1

        if namen[i] != 0 :
            for j in range(2):
                cell = int(centerg[i][j][2])
                if cell != -1:
                    center = centerg[i][j]
                    Rx = center[0] - vertexn[i][0]
                    Ry = center[1] - vertexn[i][1]
                    I_xx[i] += (Rx * Rx)
                    I_yy[i] += (Ry * Ry)
                    I_xy[i] += (Rx * Ry)
                    R_x[i] += Rx
                    R_y[i] += Ry
                    number[i] = number[i] + 1
        
        if  namen[i] == 10 and haloidn[i][-1] > 0 :
            for j in range(haloidn[i][-1]):
                cell = haloidn[i][j]
                center = centerh[cell]
                Rx = center[0] - vertexn[i][0]
                Ry = center[1] - vertexn[i][1]
                I_xx[i] += (Rx * Rx)
                I_yy[i] += (Ry * Ry)
                I_xy[i] += (Rx * Ry)
                R_x[i] += Rx
                R_y[i] += Ry
                number[i] = number[i] + 1

        D = I_xx[i]*I_yy[i] - I_xy[i]*I_xy[i]
        lambda_x[i] = (I_xy[i]*R_y[i] - I_yy[i]*R_x[i]) / D
        lambda_y[i] = (I_xy[i]*R_x[i] - I_xx[i]*R_y[i]) / D
        
    return R_x, R_y, lambda_x,lambda_y, number
