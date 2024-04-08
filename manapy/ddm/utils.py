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


def create_mpi_graph(neighbors):
    topo = COMM.Create_dist_graph_adjacent(neighbors, neighbors,
                                           sourceweights=None, destweights=None)
    return topo
    

@njit(fastmath=True)
def add(sol_1, sol_2):
    sol_1.h   += sol_2.h
    sol_1.hu  += sol_2.hu
    sol_1.hv  += sol_2.hv
    sol_1.hP  += sol_2.hP
    sol_1.hB1 += sol_2.hB1
    sol_1.hB2 += sol_2.hB2
    sol_1.Ch  += sol_2.Ch
    sol_1.Z   += sol_2.Z

    return sol_1

@njit(fastmath=True)
def minus(sol_1, sol_2):
    sol_1.h   -= sol_2.h
    sol_1.hu  -= sol_2.hu
    sol_1.hv  -= sol_2.hv
    sol_1.hP  -= sol_2.hP
    sol_1.hB1 -= sol_2.hB1
    sol_1.hB2 -= sol_2.hB2
    sol_1.Ch  -= sol_2.Ch
    sol_1.Z   -= sol_2.Z

    return sol_1

@njit(fastmath=True)
def sgn(l, eps):

    if l < -eps :
            s  = -1.
            pi = -1/l
    elif l > eps :
            s  = 1
            pi = 1/l

    else :
            s  = 0.0
            pi = 0.0  
    return s, pi    

@jit(fastmath=True)#nopython = False)
def matmul(rmatrix, matrix1, matrix2):
    
    lenmatrix1 = len(matrix1)
    lenmatrix2 = len(matrix2)
    
    for i in range(lenmatrix1):
        for j in range(len(matrix2[0])):
            for k in range(lenmatrix2):
                rmatrix[i][j] += matrix1[i][k] * matrix2[k][j]
    return rmatrix

@njit(fastmath=True)
def variables(centerc, cellidn, haloidn, vertexn, namen, centerh, centergn, halocentergn, shift, periodicn):
    
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
            number[i] += 1


        if centergn[i][0][2] != -1: 
            for j in range(len(centergn[i])):
                cell = int(centergn[i][j][-1])
                if cell != -1:
                    center = centergn[i][j]
                    Rx = center[0] - vertexn[i][0]
                    Ry = center[1] - vertexn[i][1]
                    I_xx[i] += (Rx * Rx)
                    I_yy[i] += (Ry * Ry)
                    I_xy[i] += (Rx * Ry)
                    R_x[i] += Rx
                    R_y[i] += Ry
                    number[i] = number[i] + 1


        if vertexn[i][3] == 55 or vertexn[i][3] == 56 :
            for j in range(len(periodicn[i])):
                cell = periodicn[i][j]
                if cell != -1:
                    center[0] = centerc[cell][0] + shift[cell][0]
                    center[1] = centerc[cell][1]
                    
                    Rx = center[0] - vertexn[i][0]
                    Ry = center[1] - vertexn[i][1]
                    I_xx[i] += (Rx * Rx)
                    I_yy[i] += (Ry * Ry)
                    I_xy[i] += (Rx * Ry)
                    R_x[i] += Rx
                    R_y[i] += Ry
                    number[i] += 1
      

        elif vertexn[i][3] == 57 or vertexn[i][3] == 58:
            for j in range(len(periodicn[i])):
                cell = periodicn[i][j]
                if cell != -1:
                    center[0] = centerc[cell][0]
                    center[1] = centerc[cell][1] + shift[cell][1]
                    
                    Rx = center[0] - vertexn[i][0]
                    Ry = center[1] - vertexn[i][1]
                    I_xx[i] += (Rx * Rx)
                    I_yy[i] += (Ry * Ry)
                    I_xy[i] += (Rx * Ry)
                    R_x[i] += Rx
                    R_y[i] += Ry
                    number[i] += 1
                    
        if  namen[i] == 10 :
            for j in range(len(halocentergn[i])):
                cell = int(halocentergn[i][j][-1])
                if cell != -1:
                    center[:] = halocentergn[i][j][:]
                    Rx = center[0] - vertexn[i][0]
                    Ry = center[1] - vertexn[i][1]
                    
                    I_xx[i] += (Rx * Rx)
                    I_yy[i] += (Ry * Ry)
                    I_xy[i] += (Rx * Ry)
                    R_x[i] += Rx
                    R_y[i] += Ry
                    number[i] = number[i] + 1
            
            if haloidn[i][-1] > 0:
                for j in range(haloidn[i][-1]):
                    cell = haloidn[i][j]
                    center[:] = centerh[cell][:]
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
"""
halocentergn
periodicn
"""
# TODO centertovertex
@njit(fastmath=True)
def centertovertex(w_c, wexact, w_ghost, w_halo, w_haloghost, halocentergn, centerc, centerh, cellidn, haloidn, vertexn, namen, centerg,w_node, wnode_exact, R_x, R_y, lambda_x,lambda_y, number, shift,periodicn):

    nbnode = len(vertexn)
    #center = np.zeros(3)

    for i in range(nbnode):
        for j in range(cellidn[i][-1]):
            cell = cellidn[i][j]
            center = centerc[cell]
            xdiff = center[0] - vertexn[i][0]
            ydiff = center[1] - vertexn[i][1]
            alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
            w_node[i].h  += alpha * w_c[cell].h
            w_node[i].hu += alpha * w_c[cell].hu
            w_node[i].hv += alpha * w_c[cell].hv
            w_node[i].hP += alpha * w_c[cell].hP
            w_node[i].hB1 += alpha * w_c[cell].hB1
            w_node[i].hB2 += alpha * w_c[cell].hB2
            w_node[i].Z  += alpha * w_c[cell].Z


            wnode_exact[i].h  += alpha * wexact[cell].h
            wnode_exact[i].hu += alpha * wexact[cell].hu
            wnode_exact[i].hv += alpha * wexact[cell].hv
            wnode_exact[i].hB1 += alpha * wexact[cell].hB1
            wnode_exact[i].hB2 += alpha * wexact[cell].hB2
            wnode_exact[i].Z  += alpha * wexact[cell].Z


        if centerg[i][0][2] != -1: 
            for j in range(len(centerg[i])):
                cell = int(centerg[i][j][-1])
                if cell != -1:
                    center = centerg[i][j]
                    
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]

                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                    
                    w_node[i].h   += alpha * w_ghost[cell].h
                    w_node[i].hu  += alpha * w_ghost[cell].hu
                    w_node[i].hv  += alpha * w_ghost[cell].hv
                    w_node[i].hP  += alpha * w_ghost[cell].hP
                    w_node[i].hB1 += alpha * w_ghost[cell].hB1
                    w_node[i].hB2 += alpha * w_ghost[cell].hB2
                    w_node[i].Z   += alpha * w_ghost[cell].Z

                    wnode_exact[i].h  += alpha * wexact[cellidn[i][j]].h
                    wnode_exact[i].hu += alpha * wexact[cellidn[i][j]].hu
                    wnode_exact[i].hv += alpha * wexact[cellidn[i][j]].hv
                    wnode_exact[i].hB1 += alpha * wexact[cellidn[i][j]].hB1
                    wnode_exact[i].hB2 += alpha * wexact[cellidn[i][j]].hB2
                    wnode_exact[i].Z  += alpha * wexact[cellidn[i][j]].Z
        #vertexn[i][3] == 7 or vertexn[i][3] == 8
        if vertexn[i][3] == 55 or vertexn[i][3] == 56 :
            for j in range(len(periodicn[i])):
                cell = periodicn[i][j]
                if cell != -1:
                    center = centerc[cell]
                    
                    xdiff = center[0] + shift[cell][0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                    
                    w_node[i].h  += alpha * w_c[cell].h
                    w_node[i].hu += alpha * w_c[cell].hu
                    w_node[i].hv += alpha * w_c[cell].hv
                    w_node[i].hP += alpha * w_c[cell].hP
                    w_node[i].hB1 += alpha * w_c[cell].hB1
                    w_node[i].hB2 += alpha * w_c[cell].hB2
                    w_node[i].Z  += alpha * w_c[cell].Z
                    

		    		    
	#vertexn[i][3] == 5 or vertexn[i][3] == 6	    
        elif vertexn[i][3] == 57 or vertexn[i][3] == 58:
            for j in range(len(periodicn[i])):
                cell = periodicn[i][j]
                if cell != -1:
                    center = centerc[cell]
                    
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] + shift[cell][1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
		    
                    w_node[i].h  += alpha * w_c[cell].h
                    w_node[i].hu += alpha * w_c[cell].hu
                    w_node[i].hv += alpha * w_c[cell].hv
                    w_node[i].hP += alpha * w_c[cell].hP
                    w_node[i].hB1 += alpha * w_c[cell].hB1
                    w_node[i].hB2 += alpha * w_c[cell].hB2
                    w_node[i].Z  += alpha * w_c[cell].Z
                    
				    
        if namen[i] == 10:
           
            for j in range(len(halocentergn[i])):
                cell = int(halocentergn[i][j][-1])
                if cell != -1:
                    center[:] = halocentergn[i][j][:]
                  
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                 
                    w_node[i].h   += alpha * w_haloghost[cell].h	
                    w_node[i].hu  += alpha * w_haloghost[cell].hu
                    w_node[i].hv  += alpha * w_haloghost[cell].hv
                    w_node[i].hP  += alpha * w_haloghost[cell].hP	
                    w_node[i].hB1 += alpha * w_haloghost[cell].hB1   
                    w_node[i].hB2 += alpha * w_haloghost[cell].hB2 
                    w_node[i].Z   += alpha * w_haloghost[cell].Z
                    
            if haloidn[i][-1] > 0 :
                for j in range(haloidn[i][-1]):
                    cell = haloidn[i][j]
                    center[:] = centerh[cell][:]
                  
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                                 
                    w_node[i].h   += alpha * w_halo[cell].h	
                    w_node[i].hu  += alpha * w_halo[cell].hu
                    w_node[i].hv  += alpha * w_halo[cell].hv
                    w_node[i].hP  += alpha * w_halo[cell].hP	
                    w_node[i].hB1 += alpha * w_halo[cell].hB1   
                    w_node[i].hB2 += alpha * w_halo[cell].hB2 
                    w_node[i].Z   += alpha * w_halo[cell].Z                        
                        
    return w_node, wnode_exact

@njit(fastmath=True)
def cell_gradientap(w_c, w_ghost, w_halo, w_haloghost, centerc, vertexn, nodeid, cellnid, halonid, halocenterg, namen, centerg, centerh, w_x, w_y,shift, periodicn):

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
                
                j_x = center[0]  - centerc[i][0] #+ shift[cell][0]
                j_y = center[1]  - centerc[i][1] #+ shift[cell][1]
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




@njit(fastmath=True)                     
def cell_gradient_MHD(w_c, w_ghost, w_halo, w_haloghost, centerc, cellnid, halonid, nodecid, periodicn, namen, centerg, 
                  halocenterg, vertexn, centerh, shift, nbproc, w_x, w_y):
           

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
           # j_xw += (j_x * (w_c[cell] - w_c[i] ))
           # j_yw += (j_y * (w_c[cell] - w_c[i] ))
            
        for j in range(len(periodicn[i])):
            cell = periodicn[i][j]
            if cell != -1:
                center[:] = centerc[cell][:]
                j_x = center[0] + shift[cell][0] - centerc[i][0]
                j_y = center[1] + shift[cell][1] - centerc[i][1]
                
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
               # j_xw += (j_x * (w_c[cell] - w_c[i] ))
               # j_yw += (j_y * (w_c[cell] - w_c[i] ))
                
        # if nbproc > 1:      
        for j in range(halonid[i][-1]):
            cell = halonid[i][j]
            j_x = centerh[cell][0] - centerc[i][0]
            j_y = centerh[cell][1] - centerc[i][1]
            
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
            
            #j_xw += (j_x * (w_halo[cell]  - w_c[i] ))
            #j_yw += (j_y * (w_halo[cell]  - w_c[i] ))
                
                
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





























#TODO

@njit(fastmath=True)
def barthlimiter(w_c, w_ghost, w_halo, w_x, w_y, psi, cellid, faceid, namef, halofid, centerc, centerf):
    var = "h"
    nbelement = len(w_c)
    psi[:] = 1

    for i in range(nbelement):
        w_max = w_c[var][i]
        w_min = w_c[var][i]

        for j in range(3):
            face = faceid[i][j]
            if namef[face] == 0 :#and cellid[face][1] != -1:
                w_max = max(w_max, w_c[var][cellid[face][0]], w_c[var][cellid[face][1]])
                w_min = min(w_min, w_c[var][cellid[face][0]], w_c[var][cellid[face][1]])

            elif namef[face] == 5 or namef[face] == 6 or namef[face] == 7 or namef[face] == 8:
                w_max = max(w_max, w_c[var][cellid[face][0]], w_c[var][cellid[face][1]])
                w_min = min(w_min, w_c[var][cellid[face][0]], w_c[var][cellid[face][1]])

            elif namef[face] == 1 or namef[face] == 2 or namef[face] == 3 or namef[face] == 4:
                w_max = max(w_max,  w_c[var][cellid[face][0]], w_ghost[var][face])
                w_min = min(w_min,  w_c[var][cellid[face][0]], w_ghost[var][face])
            elif namef[face] == 10:
                w_max = max(w_max,  w_c[var][cellid[face][0]], w_halo[var][halofid[face]])
                w_min = min(w_min,  w_c[var][cellid[face][0]], w_halo[var][halofid[face]])

        for j in range(3):
            face = faceid[i][j]

            r_xy = np.array([centerf[face][0] - centerc[i][0],
                            centerf[face][1] - centerc[i][1]])
            delta2 = w_x[var][i] * r_xy[0] + w_y[var][i] * r_xy[1]
            
            #TODO choice of epsilon
            if np.fabs(delta2) < 1e-10:
                psi_ij = 1.
            else:
                if delta2 > 0.:
                    value = (w_max - w_c[var][i]) / delta2
                    psi_ij = min(1., value)
                if delta2 < 0.:
                    value = (w_min - w_c[var][i]) / delta2
                    psi_ij = min(1., value)

            psi[i] = min(psi[i], psi_ij)
            
@njit(fastmath=True)
def albada(wleft, wright, w_x, w_y, center_left, center_right, lim):

    var_a = 0.
    var_b = 0.
    omega = 2./3
    epsilon = 0.
    limit = 0

    var_t = np.array([(center_right[0] - center_left[0]), (center_right[1] - center_left[1])])

    var_h = np.array([w_x.h, w_y.h])
    var_a = omega * np.dot(var_h, var_t) + (1-omega) * (wright.h - wleft.h)
    var_b = wright.h - wleft.h
    if (var_a*var_b) > 0.:
        limit = ((var_a**2 + epsilon)*var_b + (var_b**2 + epsilon)*var_a) / (var_a**2 + var_b**2)
    else:
        limit = 0.

    lim.h = limit

    var_hu = np.array([w_x.hu, w_y.hu])
    var_a = omega * np.dot(var_hu, var_t) + (1-omega) * (wright.hu - wleft.hu)
    var_b = wright.hu - wleft.hu
    if (var_a*var_b) > 0.:
        limit = ((var_a**2 + epsilon)*var_b + (var_b**2 + epsilon)*var_a) / (var_a**2 + var_b**2)
    else:
        limit = 0.

    lim.hu = limit

    var_hv = np.array([w_x.hv, w_y.hv])
    var_a = omega * np.dot(var_hv, var_t) + (1-omega) * (wright.hv - wleft.hv)
    var_b = wright.hv - wleft.hv
    if (var_a*var_b) > 0.:
        limit = ((var_a**2 + epsilon)*var_b + (var_b**2 + epsilon)*var_a) / (var_a**2 + var_b**2)
    else:
        limit = 0.

    lim.hv = limit

    var_hB1 = np.array([w_x.hB1, w_y.hB1])
    var_a = omega * np.dot(var_hB1, var_t) + (1-omega) * (wright.hB1 - wleft.hB1)
    var_b = wright.hB1 - wleft.hB1
    if (var_a*var_b) > 0.:
        limit = ((var_a**2 + epsilon)*var_b + (var_b**2 + epsilon)*var_a) / (var_a**2 + var_b**2)
    else:
        limit = 0.

    lim.hB1 = limit

    var_hB2 = np.array([w_x.hB2, w_y.hB2])
    var_a = omega * np.dot(var_hB2, var_t) + (1-omega) * (wright.hB2 - wleft.hB2)
    var_b = wright.hB2 - wleft.hB2
    if (var_a*var_b) > 0.:
        limit = ((var_a**2 + epsilon)*var_b + (var_b**2 + epsilon)*var_a) / (var_a**2 + var_b**2)
    else:
        limit = 0.

    lim.hB2 = limit

    var_Z = np.array([w_x.Z, w_y.Z])
    var_a = omega * np.dot(var_Z, var_t) + (1-omega) * (wright.Z - wleft.Z)
    var_b = wright.Z - wleft.Z
    if (var_a*var_b) > 0.:
        limit = ((var_a**2 + epsilon)*var_b + (var_b**2 + epsilon)*var_a) / (var_a**2 + var_b**2)
    else:
        limit = 0.
    lim.Z = limit
      
    
    return lim

@njit(fastmath=True)
def update(w_c, wnew, dtime, rez, src, corio, vol, cd, cr, cpsi,MGLM):
    nbelement = len(w_c)

    for i in range(nbelement):

        wnew.h[i]   = w_c.h[i]   + dtime  *  ((rez["h"][i]   + src["h"][i])/vol[i]   + corio["h"][i])
        wnew.hu[i]  = w_c.hu[i]  + dtime  *  ((rez["hu"][i]  + src["hu"][i])/vol[i]  + corio["hu"][i])
        wnew.hv[i]  = w_c.hv[i]  + dtime  *  ((rez["hv"][i]  + src["hv"][i])/vol[i]  + corio["hv"][i])
        wnew.hP[i]  = w_c.hP[i]  + dtime  *  ((rez["hP"][i]  + src["hP"][i])/vol[i]  + corio["hP"][i])
        wnew.hB1[i] = w_c.hB1[i] + dtime  *  ((rez["hB1"][i] + src["hB1"][i])/vol[i] + corio["hB1"][i])
        wnew.hB2[i] = w_c.hB2[i] + dtime  *  ((rez["hB2"][i] + src["hB2"][i])/vol[i] + corio["hB2"][i])
        wnew.Z[i]   = w_c.Z[i]   + dtime  *  ((rez["Z"][i]   + src["Z"][i])/vol[i]   + corio["Z"][i])
        wnew.Ch[i]  = rez["Ch"][i]
        
        if MGLM=="on":
             ch = cpsi #w_c.Ch[i]
             cp=np.sqrt(-dtime*pow(ch,2)/np.log(cd))
             wnew.hP[i] = np.exp(-dtime*pow(ch/cp,2))*wnew.hP[i]
             
        if MGLM=="crt":
             ch = cpsi #w_c.Ch[i]
             wnew.hP[i] = np.exp(-dtime*(ch/cr))*wnew.hP[i]        
        

    return wnew




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
                lam2 = np.fabs(u_n - B_n)
                lam3 = np.fabs(u_n + B_n)
                lam4 = np.fabs(u_n + wb)
                CPsi = max(lam1, lam2, lam3, lam4)
        c1.append(CPsi)         

    cpsi_global = max(c1)
    #print(c1)

    return cpsi_global









@njit(fastmath=True)
def time_step(w_c, cfl, normal, mesure, volume, faceid, grav, term_convective):
    nbelement =  len(faceid)
    dt_c = 1e6
    kl=1

    for i in range(nbelement):
        
        if w_c[i].h <= 0.0 :
             print(w_c[i].h)
        velson = np.sqrt(grav*w_c[i].h)
        lam = 0
        if term_convective == "non" :
            for j in range(3):
                norm = normal[faceid[i][j]]
                mesn = mesure[faceid[i][j]]
                u_n = w_c.hu[i]/w_c.h[i]*norm[0] + w_c.hv[i]/w_c.h[i]*norm[1]
                u_n = u_n/mesn
                B_n = w_c.hB1[i]/w_c.h[i]*norm[0] + w_c.hB2[i]/w_c.h[i]*norm[1]
                B_n = B_n/mesn
                wb  = np.sqrt(B_n**2 + grav*w_c[i].h)
                lam1 = np.fabs(u_n - wb)
                lam2 = np.fabs(u_n - B_n)
                lam3 = np.fabs(u_n + B_n)
                lam4 = np.fabs(u_n + wb)
                lam_convect = max(lam1, lam2, lam3, lam4)
                lam += lam_convect * mesn
                
        if kl == 1:
            for j in range(3):
                norm = normal[faceid[i][j]]
                mesn = mesure[faceid[i][j]]
                u_n = np.fabs(w_c.hu[i]/w_c.h[i]*norm[0] + w_c.hv[i]/w_c.h[i]*norm[1])
                lam_convect = u_n/mesn + velson
                lam += lam_convect * mesn

        if lam != 0:                 
            dt_c = min(dt_c, cfl * volume[i]/lam)

    dtime = np.asarray(dt_c)#np.min(dt_c))

    return dtime

@njit(fastmath=True)
def exact_solution(wexact, center, time, grav, choix):
    nbelement = len(center)
    
    choix = 90

    for i in range(nbelement):
        xcent = center[i][0]
        ycent = center[i][1]
        umax = 0.2  ;   Bmax = 0.1  ;   hmax = 1.0
        rcent = np.sqrt((xcent-time)**2 + (ycent-time)**2) 
        ee    = np.exp(1-rcent**2)
        e1    = np.exp(0.5*(1-rcent**2))
        hin   = hmax - (1./(2.0*grav))*(umax**2 - Bmax**2)*ee
        uin   = 1.0 - umax*e1*(ycent-time)
        vin   = 1.0 + umax*e1*(xcent-time)
        B1in  = -Bmax*e1*(ycent-time)
        B2in  = Bmax*e1*(xcent-time)
        wexact[i].h = hin 
        wexact[i].hu = uin     
        wexact[i].hv = vin
        wexact[i].hB1 = B1in
        wexact[i].hB2 = B2in             
        if choix == 9:
            
            c1 = 0.04
            c2 = 0.02
            theta = np.pi/6
            x0 = -20
            y0 = -10
            M = .5

            f =  -c2*((xcent - x0 - M*time*np.cos(theta))**2 + (ycent - y0 - M*time*np.sin(theta))**2)
          
            wexact[i].h = 1 - ((c1**2)/(4*grav*c2))*np.exp(2*f)
            wexact[i].hu = M * np.cos(theta) + c1*(ycent - y0 - M*time*np.sin(theta))*np.exp(f)
            wexact[i].hv = M * np.sin(theta) - c1*(xcent - x0 - M*time*np.cos(theta))*np.exp(f)

            wexact[i].hu = wexact[i].hu * wexact[i].h
            wexact[i].hv = wexact[i].hv * wexact[i].h
    return wexact

#from numba.typed import Dict

#TODO initialisation
@njit(fastmath=True)
def initialisation(w_c, center, choix):

    nbelements = len(center)




    if choix == 1:  # Orszag–Tang-like turbulence problem
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            ha = 1./4
            ua = 1 + 0.5*np.sin(np.pi*ycent) + ha*np.cos(np.pi*xcent)
            va = 1 + ha*np.sin(np.pi*xcent) + 0.5*np.cos(np.pi*ycent)

            
            
            w_c.h[i]   = ha
            w_c.hu[i]  = ha*ua
            w_c.hv[i]  = ha*va
            w_c.hB1[i] = ha*(0.5)
            w_c.hB2[i] = ha*ua  

            w_c.hP[i] = 0.   
            w_c.Ch[i] = 0.







  

    if choix == 2:  # One dimensional dam break problem in SWMHD
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            if xcent <= 0 :     
                
                w_c.hu[i] = 0.
                w_c.h[i]  = 1.
                w_c.hv[i] = 0.
                w_c.Z[i]  = 0.
                w_c.hB1[i] = 1.
                w_c.hB2[i] = 0.
            else:
                
                w_c.Z[i]  = 0. 
                w_c.h[i]  = 2.  
                w_c.hu[i] = 0.
                w_c.hv[i] = 0.
                w_c.hB1[i] = 1.
                w_c.hB2[i] = 2.
            w_c.hP[i] = 0.
            w_c.Ch[i] = 0.
            
            

    elif choix == 3: #Two dimensional explosion problem
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            if np.sqrt(xcent**2 + ycent**2) <= 0.3 :     
                w_c.h[i]   = 1.
                w_c.hu[i]  = 0.
                w_c.hv[i]  = 0.
                w_c.Z[i]   = 0.
                w_c.hB1[i] = 0.1
                w_c.hB2[i] = 0.0
            else:
                w_c.h[i]   = 0.1
                w_c.Z[i]   = 0.   
                w_c.hu[i]  = 0.
                w_c.hv[i]  = 0.
                w_c.hB1[i] = 0.1
                w_c.hB2[i] = 0.0

            w_c.hP[i] = 0.   
            w_c.Ch[i] = 0.

            
    elif choix == 4: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]

            g = 1.0  ;   umax = 0.2  ;   Bmax = 0.1  ;   hmax = 1.0
            
            rcent = np.sqrt(xcent**2 + ycent**2)            
            ee    = np.exp(1-rcent**2)
            e1    = np.exp(0.5*(1-rcent**2))
            hin   = hmax - 0.5*(umax**2 - Bmax**2)*ee
            uin   = 1.0 - umax*e1*ycent
            vin   = 1.0 + umax*e1*xcent
            B1in  = -Bmax*e1*ycent
            B2in  = Bmax*e1*xcent
            
            w_c.h[i]  = hin
            w_c.hu[i] = hin*uin
            w_c.hv[i] = hin*vin
            w_c.hB1[i] = hin*B1in
            w_c.hB2[i] = hin*B2in

            w_c.hP[i] = 0.   
            w_c.Ch[i] = 0. 
 
         

    elif choix == 5: #h gaussienne (c-propriété)

        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            

            w_c.Z[i] = 0.8 * np.exp(-5*(xcent - 1)**2 - 50* (ycent - 0.5)**2) 


            w_c.h[i]   = 1-w_c.Z[i] 
            w_c.hu[i]  = 0.
            w_c.hv[i]  = 0.
            w_c.hB1[i] = 0.
            w_c.hB2[i] = 0.
            w_c.hP[i]  = 0.
            w_c.Ch[i] = 0. 



    elif choix == 6: #h gaussienne (c-propriété)

        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            

            w_c.Z[i] = 0.8 * np.exp(-5*(xcent - 1)**2 - 50* (ycent - 0.5)**2) 

            
            if 0.05 < xcent and xcent < 0.5 :
            	   hint = 1. - w_c.Z[i] #+0.01
            else :
            	   hint = 1 - w_c.Z[i]

            w_c.h[i] = hint
            w_c.hu[i] = 0.
            w_c.hv[i] = 0.
            w_c.hB1[i] = 1
            w_c.hB2[i] = 1e-4
            w_c.hP[i] = 0.




    elif choix == 8:  # Orszag–Tang-like turbulence problem
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            gamma = 1.667
            hint = 25/9# pow(gamma,2)
            
            
            w_c.h[i]   =  hint
            w_c.hu[i]  = -hint*np.sin(ycent) 
            w_c.hv[i]  =  hint*np.sin(xcent) 
            w_c.hB1[i] = -hint*np.sin(ycent) 
            w_c.hB2[i] =  hint*np.sin(2*xcent)   

            w_c.hP[i] = 0.   
            w_c.Ch[i] = 0.


    return w_c



#TODO ghost
@njit
def ghost_value(w_c, w_ghost, cellid, name, normal, mesure, time, center, centerf,cellnid, w_x, w_y):


    nbface = len(cellid)
  
        
    for i in range(nbface):
        w_ghost[i] = w_c[cellid[i][0]]
        #xcent = center[i][0]
        #xcent = center[i][1]
        xc    = center[cellid[i][0]][0]
        yc    = center[cellid[i][0]][1]
        xe    = centerf[i][0]
        ye    = centerf[i][1]

        if ( name[i] == 41 or name[i] == 42 or name[i] == 43 or name[i] == 44 ):         
          
            w_ghost[i].hu = w_c[cellid[i][0]].hu
            w_ghost[i].hv = w_c[cellid[i][0]].hv
            xcent = normal[i][0]
            ycent = normal[i][1]
            
            umax = 0.2  ;   Bmax = 0.1  ;   hmax = 1.0
            rcent = np.sqrt((xcent-time)**2 + (ycent-time)**2)
            ee    = np.exp(1-rcent**2) 
            e1    = np.exp(0.5*(1-rcent**2))
            hin   = hmax - (1./(2.0))*(umax**2 - Bmax**2)*ee
            uin   = 1.0 - umax*e1*(ycent-time)
            vin   = 1.0 + umax*e1*(xcent-time)
            B1in  = -Bmax*e1*(ycent-time)
            B2in  = Bmax*e1*(xcent-time)
            w_ghost[i].h = hin 
            w_ghost[i].hu = hin*uin
            w_ghost[i].hv = hin*vin
            w_ghost[i].hB1 = hin*B1in
            w_ghost[i].hB2 = hin*B2in
            w_ghost[i].hP = w_c[cellid[i][0]].hP
     
        if ( name[i] == 1 or name[i] == 2 ):         
          
            w_ghost[i].hu = w_c[cellid[i][0]].hu
            w_ghost[i].hv = w_c[cellid[i][0]].hv
            w_ghost[i].hP = w_c[cellid[i][0]].hP
            w_ghost[i].h  = w_c[cellid[i][0]].h 
            w_ghost[i].Z  = w_c[cellid[i][0]].Z
            w_ghost[i].hB1 = w_c[cellid[i][0]].hB1 #+ (xe - xc)*w_x[cellid[i][0]].hB1 + (ye - yc)*w_y[cellid[i][0]].hB1
            w_ghost[i].hB2 = w_c[cellid[i][0]].hB2 #+ (xe - xc)*w_x[cellid[i][0]].hB2 + (ye - yc)*w_y[cellid[i][0]].hB2
            
        if (name[i] == 3 or name[i] == 4):
            
            u_i = w_c[cellid[i][0]].hu/w_c[cellid[i][0]].h
            v_i = w_c[cellid[i][0]].hv/w_c[cellid[i][0]].h
            B1_i = w_c[cellid[i][0]].hB1/w_c[cellid[i][0]].h
            B2_i = w_c[cellid[i][0]].hB2/w_c[cellid[i][0]].h
            
 
            s_n = normal[i] / mesure[i]
            
            B1n_g = B1_i*(s_n[1]*s_n[1] - s_n[0]*s_n[0]) - 2.0*B2_i*s_n[0]*s_n[1]
            B2n_g = B2_i*(s_n[0]*s_n[0] - s_n[1]*s_n[1]) - 2.0*B1_i*s_n[0]*s_n[1]
        
            u_g = u_i*(s_n[1]*s_n[1] - s_n[0]*s_n[0]) - 2.0*v_i*s_n[0]*s_n[1]
            v_g = v_i*(s_n[0]*s_n[0] - s_n[1]*s_n[1]) - 2.0*u_i*s_n[0]*s_n[1]
            
         # face.cellid[i][0] au bord c'est la cell gauche [i][1] à droite
            # cells.cellf_id[face.cellid[i][0]][0:3]
            w_ghost[i].h = w_c[cellid[i][0]].h
            w_ghost[i].Z = w_c[cellid[i][0]].Z
            w_ghost[i].hP = w_c[cellid[i][0]].hP
            w_ghost[i].hB1 = w_c[cellid[i][0]].h * B1n_g
            w_ghost[i].hB2 = w_c[cellid[i][0]].h * B2n_g
            w_ghost[i].hu = w_c[cellid[i][0]].h * u_g
            w_ghost[i].hv = w_c[cellid[i][0]].h * v_g
            


    return w_ghost


#TODO ghost
@njit
def ghost_valueff(w_c, w_ghost, cellid, name, normal, mesure, time, center, centerf,cellnid, w_x, w_y):


    nbface = len(cellid)
  
        
    for i in range(nbface):

        theta = np.pi/6
        M = .5
        c1 = 0.04
        c2 = 0.02
        grav = 1.
        x0 = -20
        y0 = -10

        w_ghost[i] = w_c[cellid[i][0]]

        if  name[i] == 1:   
            
            xcent = normal[i][0]
            ycent = normal[i][1]
            
            f =  -c2*((-50 - x0 - M*time*np.cos(theta))**2 + (ycent - y0 - M*time*np.sin(theta))**2)
            
            
            w_ghost[i].h = 1 - ((c1**2)/(4*grav*c2))*np.exp(2*f)
          
            u = M * np.cos(theta) + c1*(ycent - y0 - M*time*np.sin(theta))*np.exp(f)
            v = M * np.sin(theta) - c1*(-50 - x0 - M*time*np.cos(theta))*np.exp(f)
            
            w_ghost[i].hu = w_ghost[i].h * u
             
            w_ghost[i].hv = w_ghost[i].h * v

        if  name[i] == 2 :
             xcent = normal[i][0]
             ycent = normal[i][1]
            
             f =  -c2*((50 - x0 - M*time*np.cos(theta))**2 + (ycent - y0 - M*time*np.sin(theta))**2)
            
             w_ghost[i].h = 1 - ((c1**2)/(4*grav*c2))*np.exp(2*f)
          
             u = M * np.cos(theta) + c1*(ycent - y0 - M*time*np.sin(theta))*np.exp(f)
             v = M * np.sin(theta) - c1*(50 - x0 - M*time*np.cos(theta))*np.exp(f)
            
             w_ghost[i].hu = w_ghost[i].h * u
             w_ghost[i].hv = w_ghost[i].h * v

        if  name[i] == 3 :
            
             xcent = normal[i][0] 
             ycent = normal[i][1]
            
             f =  -c2*((xcent - x0 - M*time*np.cos(theta))**2 + (50 - y0 - M*time*np.sin(theta))**2)
            
             w_ghost[i].h = 1 - ((c1**2)/(4*grav*c2))*np.exp(2*f)
          
             u = M * np.cos(theta) + c1*(50 - y0 - M*time*np.sin(theta))*np.exp(f)
             v = M * np.sin(theta) - c1*(xcent - x0 - M*time*np.cos(theta))*np.exp(f)
            
             w_ghost[i].hu = w_ghost[i].h * u
             w_ghost[i].hv = w_ghost[i].h * v 

        if  name[i] == 4 :

               xcent = normal[i][0] 
               ycent = normal[i][1]
             
               f =  -c2*((xcent - x0 - M*time*np.cos(theta))**2 + (-50 - y0 - M*time*np.sin(theta))**2)
             
               w_ghost[i].h = 1 - ((c1**2)/(4*grav*c2))*np.exp(2*f)
           
               u = M * np.cos(theta) + c1*(-50 - y0 - M*time*np.sin(theta))*np.exp(f)
               v = M * np.sin(theta) - c1*(xcent - x0 - M*time*np.cos(theta))*np.exp(f)
             
               w_ghost[i].hu = w_ghost[i].h * u
               w_ghost[i].hv = w_ghost[i].h * v
        
    return w_ghost



@njit
def update_haloghost_value_MHD(w_halo, w_haloghost, haloghostcenter, vertexn, namen, time,):
    
    nbnodes = len(vertexn)
    
    for i in range(nbnodes):
        if namen[i] == 10:
            for j in range(len(haloghostcenter[i])):
                
                if haloghostcenter[i][j][-1] != -1:
                    cellhalo  = int(haloghostcenter[i][j][-3])
                    cellghost = int(haloghostcenter[i][j][-1])
                    
                    w_haloghost[cellghost].h   = w_halo[cellhalo].h
                    w_haloghost[cellghost].hu  = w_halo[cellhalo].hu
                    w_haloghost[cellghost].hv  = w_halo[cellhalo].hv
                    w_haloghost[cellghost].hB1  = w_halo[cellhalo].hB1
                    w_haloghost[cellghost].hB2  = w_halo[cellhalo].hB2
                    w_haloghost[cellghost].Z   = w_halo[cellhalo].Z
    
    return w_haloghost
                        

def norml1(w_c, vol, cells, nodes, w_x, w_y, ET):


    cells  = np.array(cells)
    nbelement = len(w_c)
    divhB = np.zeros(len(w_c))
    ETL   = np.zeros(len(w_c))
    eng=0.0


    for i in range(nbelement):

        divhB[i] = (w_x["hB1"][i] + w_y["hB2"][i])*vol[i]
        #ETL[i]   = ET[i]*vol[i]
        #eng += np.fabs(ET[i])*vol[i]
        


    dErrorL2 = np.linalg.norm(divhB,ord=1)
    div = dErrorL2
    #eng = np.linalg.norm(ETL,ord=2)
    
    return div#, eng

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
        


    dErrorL2 = np.linalg.norm(divhB,ord=2)
    div = dErrorL2
    #eng = np.linalg.norm(ETL,ord=2)
    
    return div, np.sqrt(eng)


def normlinf(w_c, vol, cells, nodes, w_x, w_y, ET):


    cells  = np.array(cells)
    nbelement = len(w_c)
    divhB = np.zeros(len(w_c))
    ETL   = np.zeros(len(w_c))
    eng=0.0


    for i in range(nbelement):

        divhB[i] = (w_x["hB1"][i] + 0.0*w_y["hB2"][i])*vol[i]
        #ETL[i]   = ET[i]*vol[i]     


    dErrorL2 = np.linalg.norm(divhB,ord=np.inf)
    div = dErrorL2
    #eng = np.linalg.norm(ETL,ord=2)
    
    return div#, max(ETL)



def Total_Energy(w_c, grav, correction):

    nbelement = len(w_c)
    Et = np.zeros(nbelement) # Total_Energy
    Ec = np.zeros(nbelement) # Kinetic Energy
    Ep = np.zeros(nbelement) # Potential Energy
    Em = np.zeros(nbelement) # Magnetic energy
    EN = 0.0
    
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
        if correction == "on" :
            Et[i] = Ec[i] + Em[i] + Ep[i] + 0.5*pow(w_c["hP"][i] , 2)
        else :
            Et[i] = Ec[i] + Em[i] + Ep[i]  
    return Et

def Total_En(w_c, grav, correction, vol):

    nbelement = len(w_c)
    Et = np.zeros(nbelement) # Total_Energy
    Ec = np.zeros(nbelement) # Kinetic Energy
    Ep = np.zeros(nbelement) # Potential Energy
    Em = np.zeros(nbelement) # Magnetic energy
    EN = 0.0
    
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
        if correction == "on" :
            Et[i] = Ec[i] + Em[i] + Ep[i] + 0.5*pow(w_c["hP"][i] , 2)
        else :
            Et[i] = Ec[i] + Em[i] + Ep[i]  
            EN += Et[i]*vol[i]
    return EN










def save_paraview_results(w_c, wexact, niter, miter, time, dtime, rank, size, vol, cells, nodes, w_x, w_y, tf, ET):


    elements = {"triangle": cells}
    points = []
    for i in nodes:
        points.append([i[0], i[1], i[2]])
    
    cells  = np.array(cells)
    points = np.array(points)
    nbelement = len(w_c)
    Error_h = np.zeros(len(w_c))
    Error_u = np.zeros(len(w_c))
    Error_v = np.zeros(len(w_c))
    Error_B1 = np.zeros(len(w_c))
    Error_B2 = np.zeros(len(w_c))
    he = np.zeros(len(w_c))
    ue = np.zeros(len(w_c))
    divhB = np.zeros(len(w_c))
    
#    divhBd = np.zeros(len(w_c))
#    TimeS = []
#    L2div = []

    for i in range(nbelement):

        divhB[i] =(w_x["hB1"][i] + w_y["hB2"][i])*vol[i] #DIV[i]
        Error_h[i] = np.fabs(w_c["h"][i] - wexact["h"][i])*vol[i]
        Error_u[i] = np.fabs(w_c["hu"][i]/w_c["h"][i]   - wexact["hu"][i])*vol[i]
        Error_v[i] = np.fabs(w_c["hv"][i]/w_c["h"][i]   - wexact["hv"][i])*vol[i]
        Error_B1[i] = np.fabs(w_c["hB1"][i]/w_c["h"][i] - wexact["hB1"][i])*vol[i]
        Error_B2[i] = np.fabs(w_c["hB2"][i]/w_c["h"][i] - wexact["hB2"][i])*vol[i]
#        divhBd[i] =(w_x["hB1"][i] + w_y["hB2"][i])#*vol[i]
#        Error_h[i] = np.fabs(w_c["h"][i] - wexact["h"][i])*vol[i]
#        Error_u[i] = np.fabs(w_c["hu"][i]/w_c["h"][i] - wexact["hu"][i]/wexact["h"][i])*vol[i]
#        he[i]=np.fabs(wexact["h"][i])*vol[i]
#        ue[i]=np.fabs(wexact["hu"][i]/wexact["h"][i])*vol[i]     

#    ErrorL1 =np.linalg.norm(divhB,ord=1)
    dErrorL2 = np.linalg.norm(divhB,ord=2)
    div = dErrorL2
#    ErrorLinf =np.linalg.norm(divhB,ord=np.inf)

    #ErrorL1 =np.linalg.norm(w_c["hP"]/w_c["h"],ord=1)
    ErrorL2 =np.linalg.norm(w_c["hP"]/w_c["h"],ord=2)
    #ErrorLinf =np.linalg.norm(w_c["hP"]/w_c["h"],ord=np.inf)
    
    ErrorL1_h  = np.linalg.norm(Error_h,ord=1)/np.linalg.norm(wexact["h"],ord=1)
    ErrorL1_u  = np.linalg.norm(Error_u,ord=1)/np.linalg.norm(wexact["hu"],ord=1)
    ErrorL1_v  = np.linalg.norm(Error_v,ord=1)/np.linalg.norm(wexact["hv"],ord=1)
    ErrorL1_B1 = np.linalg.norm(Error_B1,ord=1)/np.linalg.norm(wexact["hB1"],ord=1)
    ErrorL1_B2 = np.linalg.norm(Error_B2,ord=1)/np.linalg.norm(wexact["hB2"],ord=1)
    
    data = {"h" : w_c["h"], "u" : w_c["hu"]/w_c["h"], "v": w_c["hv"]/w_c["h"], "B1": w_c["hB1"]/w_c["h"], "B2": w_c["hB2"]/w_c["h"], "Z": w_c["Z"], "h+Z": w_c["h"] + w_c["Z"], "uex" : wexact["hu"], "divhB": divhB, "B": np.sqrt((w_c["hB1"]/w_c["h"])**2 + (w_c["hB2"]/w_c["h"])**2), "vex": wexact["hv"]}
   
    if len(w_c) == len(cells):
        data = {"h": data, "u":data, "v": data, "B1": data, "B2": data, "Z":data, "h+Z":data,  "uex":data, "divhB":data, "ET":data, "vex":data}

    maxhZZ = max(w_c["h"]+w_c["Z"])
    minhZZ = min(w_c["h"]+w_c["Z"])
    maxhZ  = max(np.fabs(divhB))
    minhZ  = min(np.fabs(divhB))
    maxu   = max(w_c["hu"]/w_c["h"])
#    maxhZ1 = max(np.fabs(divhBd))
#    minhZ1 = min(np.fabs(divhBd))    
    integral_maxhZ = np.zeros(1)
    integral_minhZ = np.zeros(1)
    integral_maxhZZ = np.zeros(1)
    integral_minhZZ = np.zeros(1)
#    integral_maxhZ1 = np.zeros(1)
#    integral_minhZ1 = np.zeros(1)
    integral_maxu = np.zeros(1, dtype=np.double)

    COMM.Reduce(maxhZ, integral_maxhZ, MPI.MAX, 0)
    COMM.Reduce(minhZ, integral_minhZ, MPI.MIN, 0)
    COMM.Reduce(maxhZZ, integral_maxhZZ, MPI.MAX, 0)
    COMM.Reduce(minhZZ, integral_minhZZ, MPI.MIN, 0)
#    COMM.Reduce(maxhZ1, integral_maxhZ1, MPI.MAX, 0)
#    COMM.Reduce(minhZ1, integral_minhZ1, MPI.MIN, 0)    
    COMM.Reduce(maxu, integral_maxu, MPI.MAX, 0)
    if rank == 0:
        print(" **************************** Computing ****************************")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Saving Results $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print("Iteration = ", niter, "time = ", np.float32(time), "time step = ", np.float32(dtime))
        #print("max h+Z =", np.float32(integral_maxhZ[0]), "max u =", np.float32(integral_maxu[0]))
        #print("max divhB =", np.float32(integral_maxhZ[0]), "max u =", np.float32(integral_maxu[0]))
        #print("min divhB =", np.float32(integral_minhZ[0]))
        print("max h+Z =", np.float32(integral_maxhZZ[0]), "min h+Z =", np.float32(integral_minhZZ[0]))
        #print("Error_L1 for h =", ErrorL1_h)
        #print("Error_L1 for h =", ErrorL1_u, "Error_L1 for v =", ErrorL1_v)
        #print("Error_L1 for B1 =", ErrorL1_B1, "Error_L1 for B2 =", ErrorL1_B2)

#        print("max divhBd =", np.float32(integral_maxhZ1[0]))
#        print("min divhBd =", np.float32(integral_minhZ1[0]))

        #print("Error_L1 for divhB =", ErrorL1, "Error_L2 for divhB =", ErrorL2,"Error_Linf for divhB =", ErrorLinf)
        #print("Error_L1 for Psi =", ErrorL1, "Error_L2 for Psi =", ErrorL2,"Error_Linf for Psi =", ErrorLinf)
        #print("Error_L2 for divhB =", dErrorL2, "Error_L2 for Psi =", ErrorL2)
#        TimeS.append(time)
#        L2div.append(dErrorL2)
        div = dErrorL2
#        if (time > tf-0.001 and time < tf + 0.1) :
#             #open('Temps.txt','wb')
#             print("temps=", TimeS)
#             print("div=", L2div)
#             with open("Temps.txt", "w") as f :
#                  np.savetxt(f, TimeS, fmt='%.4f')
#             with open("divhb.txt", "w") as f :
#                  np.savetxt(f, L2div, fmt='%.4f')
             #open('divhb.txt','wb')
             #np.savetxt(divhb.txt, L2div, fmt='%.2f')

    
    if len(w_c) == len(cells):
        meshio.write_points_cells("results/visu"+str(rank)+"-"+str(miter)+".vtu",
                                  points, elements, cell_data=data, file_format="vtu")
    else:
        meshio.write_points_cells("results/visu"+str(rank)+"-"+str(miter)+".vtu",
                                  points, elements, point_data=data, file_format="vtu")

    if(rank == 0 and size > 1):
        with open("results/visu"+str(miter)+".pvtu", "a") as text_file:
            text_file.write("<?xml version=\"1.0\"?>\n")
            text_file.write("<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
            text_file.write("<PUnstructuredGrid GhostLevel=\"0\">\n")
            text_file.write("<PPoints>\n")
            text_file.write("<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"binary\"/>\n")
            text_file.write("</PPoints>\n")
            text_file.write("<PCells>\n")
            text_file.write("<PDataArray type=\"Int64\" Name=\"connectivity\" format=\"binary\"/>\n")
            text_file.write("<PDataArray type=\"Int64\" Name=\"offsets\" format=\"binary\"/>\n")
            text_file.write("<PDataArray type=\"Int64\" Name=\"types\" format=\"binary\"/>\n")
            text_file.write("</PCells>\n")
            if len(w_c) == len(cells):
                text_file.write("<PCellData Scalars=\"h\">\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"h\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"u\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"v\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"B1\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"B2\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"h+Z\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"divhB\" format=\"binary\"/>\n")
                text_file.write("</PCellData>\n")
            else:
                text_file.write("<PPointData Scalars=\"h\">\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"h\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"u\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"v\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"B1\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"B2\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"h+Z\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"divhB\" format=\"binary\"/>\n")
                text_file.write("</PPointData>\n")
            for i in range(size):
                name1 = "visu"
                bu1 = [10]
                bu1 = str(i)
                name1 += bu1
                name1 += "-"+str(miter)
                name1 += ".vtu"
                text_file.write("<Piece Source=\""+str(name1)+"\"/>\n")
            text_file.write("</PUnstructuredGrid>\n")
            text_file.write("</VTKFile>")
