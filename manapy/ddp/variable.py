#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 19:04:50 2020

@author: kissami
"""
import meshio
from mpi4py import MPI
from numpy import (array, arange, zeros, sin, cos, pi, fabs, sqrt, double, float16, 
                   float64, exp, linspace, ones, timedelta64, datetime64, zeros)
from numba import njit

COMM = MPI.COMM_WORLD

     



@njit(fastmath=True)    
def variables_MHD(centerc:'float[:,:]', cellidn:'int[:,:]', haloidn:'int[:,:]', periodicn:'int[:,:]', 
              vertexn:'float[:,:]', namen:'int[:]', 
              centergn:'float[:,:,:]', halocentergn:'float[:,:,:]', centerh:'float[:,:]', nbproc:'int',  
              R_x:'float[:]', R_y:'float[:]', lambda_x:'float[:]', 
              lambda_y:'float[:]', number:'int[:]', shift:'float[:,:]'):
    

      nbnode = len(R_x)
        
      I_xx = zeros(nbnode)
      I_yy = zeros(nbnode)
      I_xy = zeros(nbnode)
      center = zeros(3)
   
      for i in range(nbnode):
        for j in range(cellidn[i][-1]):
            center[:] = centerc[cellidn[i][j]][:]
            Rx = center[0] - vertexn[i][0]
            Ry = center[1] - vertexn[i][1]
            I_xx[i] += (Rx * Rx)
            I_yy[i] += (Ry * Ry)
            I_xy[i] += (Rx * Ry)
            R_x[i] += Rx
            R_y[i] += Ry
            number[i] += 1
            
        #ghost boundary (old vertex names)
        # if vertexn[i][3] == 1 or vertexn[i][3] == 2 or vertexn[i][3] == 3 or vertexn[i][3] == 4:
        if centergn[i][0][2] != -1: 
            for j in range(len(centergn[i])):
                cell = int(centergn[i][j][-1])
                if cell != -1:
                    center[:] = centergn[i][j][:]
                    Rx = center[0] - vertexn[i][0]
                    Ry = center[1] - vertexn[i][1]
                    #Rx = centergn[i][j][0] - vertexn[i][0]
                    #Ry = centergn[i][j][1] - vertexn[i][1]
                    I_xx[i] += (Rx * Rx)
                    I_yy[i] += (Ry * Ry)
                    I_xy[i] += (Rx * Ry)
                    R_x[i] += Rx
                    R_y[i] += Ry
                    number[i] = number[i] + 1
        
        #periodic boundary old vertex names)
        if vertexn[i][3] == 5 or vertexn[i][3] == 6 :
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
                    
        elif vertexn[i][3] == 7 or vertexn[i][3] == 8:
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
# TODO centertovertex
@njit(fastmath=True)    
def centertovertex(h_c:'float[:]', hu_c:'float[:]', hv_c:'float[:]', hc_c:'float[:]', Z_c:'float[:]', 
                   h_ghost:'float[:]', hu_ghost:'float[:]', hv_ghost:'float[:]', hc_ghost:'float[:]', Z_ghost:'float[:]',
                   h_halo:'float[:]', hu_halo:'float[:]', hv_halo:'float[:]', hc_halo:'float[:]', Z_halo:'float[:]',
                   h_haloghost:'float[:]', hu_haloghost:'float[:]', hv_haloghost:'float[:]', hc_haloghost:'float[:]', 
                   Z_haloghost:'float[:]',
                   centerc:'float[:,:]', centerh:'float[:,:]', cellidn:'int[:,:]', haloidn:'int[:,:]', 
                   vertexn:'float[:,:]', namen:'int[:]', centergn:'float[:,:,:]', halocentergn:'float[:,:,:]', 
                   R_x:'float[:]', R_y:'float[:]', lambda_x:'float[:]',lambda_y:'float[:]', number:'float[:]', 
                   alphanbproc:'int',
                   h_n:'float[:]', hu_n:'float[:]', hv_n:'float[:]', hc_n:'float[:]', Z_n:'float[:]'):
    
    h_n[:] = 0.; hu_n[:] = 0.; hv_n[:] = 0.; hc_n[:] = 0.; Z_n[:] = 0.
    
    nbnode = len(vertexn)
    center = zeros(2)
    
    for i in range(nbnode):
        for j in range(cellidn[i][-1]):
            cell = cellidn[i][j]
            center[:] = centerc[cell][0:2]
           
            xdiff = center[0] - vertexn[i][0]
            ydiff = center[1] - vertexn[i][1]
            alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
       
            h_n[i]  += alpha * h_c[cell]
            hu_n[i] += alpha * hu_c[cell]
            hv_n[i] += alpha * hv_c[cell]
            hc_n[i] += alpha * hc_c[cell]
            Z_n[i]  += alpha * Z_c[cell]
                
        if namen[i] != 0 :
           for j in range(len(centergn[i])):
                cell = int(centergn[i][j][-1])
                if cell != -1:
                    center[:] = centergn[i][j][0:2]
                    
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                
                    h_n[i]  += alpha * h_ghost[cell]
                    hu_n[i] += alpha * hu_ghost[cell]
                    hv_n[i] += alpha * hv_ghost[cell]
                    hc_n[i] += alpha * hc_ghost[cell]
                    Z_n[i]  += alpha * Z_ghost[cell]
                
        if namen[i] == 10 and haloidn[i][-1] > 0 :
            for j in range(haloidn[i][-1]):
                cell = haloidn[i][j]
                center[:] = centerh[cell][0:2]
              
                xdiff = center[0] - vertexn[i][0]
                ydiff = center[1] - vertexn[i][1]
                alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
              
                h_n[i]  += alpha * h_halo[cell]
                hu_n[i] += alpha * hu_halo[cell]
                hv_n[i] += alpha * hv_halo[cell]
                hc_n[i] += alpha * hc_halo[cell]
                Z_n[i]  += alpha * Z_halo[cell]
                
        if namen[i] == 10 :
            for j in range(len(halocentergn[i])):
                cell = int(halocentergn[i][j][-1])
                if cell != -1:
                    center[:] = halocentergn[i][j][0:2]
                  
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                 
                    h_n[i]  += alpha * h_haloghost[cell]
                    hu_n[i] += alpha * hu_haloghost[cell]
                    hv_n[i] += alpha * hv_haloghost[cell]
                    hc_n[i] += alpha * hc_haloghost[cell]
                    Z_n[i]  += alpha * Z_haloghost[cell]
                        


