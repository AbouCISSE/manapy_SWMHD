#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 19:04:50 2020

@author: kissami
"""
import meshio
from mpi4py import MPI
import numpy as np
from numba import njit, jit

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
    sol_1.hB1 += sol_2.hB1
    sol_1.hB2 += sol_2.hB2
    sol_1.Z   += sol_2.Z

    return sol_1

@njit(fastmath=True)
def minus(sol_1, sol_2):
    sol_1.h   -= sol_2.h
    sol_1.hu  -= sol_2.hu
    sol_1.hv  -= sol_2.hv
    sol_1.hB1 -= sol_2.hB1
    sol_1.hB2 -= sol_2.hB2
    sol_1.Z   -= sol_2.Z

    return sol_1

@njit(fastmath=True)
def sgn(l, eps):
    if np.fabs(l) < eps :
            s  = 0.
            pi = 0.
    else :
            s  = l/np.fabs(l)
            pi = s/l
  
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
def variables(centerc, cellidn, haloidn, vertexn, namen, centerg, centerh):
     
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

# TODO centertovertex
@njit(fastmath=True)
def centertovertex(w_c, w_ghost, w_halo, centerc, centerh, cellidn, haloidn, vertexn, namen, centerg,
                   w_node, wexact, wnode_exact, R_x, R_y, lambda_x,lambda_y, number):

    nbnode = len(vertexn)

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
            w_node[i].hB1 += alpha * w_c[cell].hB1
            w_node[i].hB2 += alpha * w_c[cell].hB2
            w_node[i].Z  += alpha * w_c[cell].Z
            
            #Exact solution
            wnode_exact[i].h  += alpha * wexact[cell].h
            wnode_exact[i].hu += alpha * wexact[cell].hu
            wnode_exact[i].hv += alpha * wexact[cell].hv
            wnode_exact[i].hB1 += alpha * wexact[cell].hB1
            wnode_exact[i].hB2 += alpha * wexact[cell].hB2
            wnode_exact[i].Z  += alpha * wexact[cell].Z
       
        if namen[i] != 0 :
           for j in range(2):
                cell = int(centerg[i][j][2])
                if cell != -1:
                    center = centerg[i][j]
                    ghost = np.int(centerg[i][j][2])
                    xdiff = center[0] - vertexn[i][0]
                    ydiff = center[1] - vertexn[i][1]
                    alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                    w_node[i].h  += alpha * w_ghost[ghost].h
                    w_node[i].hu += alpha * w_ghost[ghost].hu
                    w_node[i].hv += alpha * w_ghost[ghost].hv
                    w_node[i].hB1 += alpha * w_ghost[ghost].hB1
                    w_node[i].hB2 += alpha * w_ghost[ghost].hB2                    
                    w_node[i].Z  += alpha * w_ghost[ghost].Z
                    
                    #Exact solution
                    wnode_exact[i].h  += alpha * wexact[cellidn[i][j]].h
                    wnode_exact[i].hu += alpha * wexact[cellidn[i][j]].hu
                    wnode_exact[i].hv += alpha * wexact[cellidn[i][j]].hv
                    wnode_exact[i].hB1 += alpha * wexact[cellidn[i][j]].hB1
                    wnode_exact[i].hB2 += alpha * wexact[cellidn[i][j]].hB2
                    wnode_exact[i].Z  += alpha * wexact[cellidn[i][j]].Z
                
        if namen[i] == 10 and haloidn[i][-1] > 0 :
            for j in range(haloidn[i][-1]):
                cell = haloidn[i][j]
                center = centerh[cell]
                xdiff = center[0] - vertexn[i][0]
                ydiff = center[1] - vertexn[i][1]
                alpha = (1. + lambda_x[i]*xdiff + lambda_y[i]*ydiff)/(number[i] + lambda_x[i]*R_x[i] + lambda_y[i]*R_y[i])
                w_node[i].h  += alpha * w_halo[cell].h
                w_node[i].hu += alpha * w_halo[cell].hu
                w_node[i].hv += alpha * w_halo[cell].hv
                w_node[i].hB1 += alpha * w_halo[cell].hB1
                w_node[i].hB2 += alpha * w_halo[cell].hB2
                w_node[i].Z  += alpha * w_halo[cell].Z
                        
                        
    return w_node, wnode_exact

@njit(fastmath=True)
def cell_gradient(w_c, w_ghost, w_halo, centerc, nodeidc, cellnid, halonid, namen, centerg, centerh,
                  w_x, w_y):

    nbelement = len(centerc)
    for i in range(nbelement):
        i_xx = 0.
        i_yy = 0.
        i_xy = 0.
        j_xh = 0.
        j_yh = 0

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
        
        if SIZE > 1:      
            for j in range(halonid[i][-1]):
                cell = halonid[i][j]
                j_x = centerh[cell][0] - centerc[i][0]
                j_y = centerh[cell][1] - centerc[i][1]
                i_xx += pow(j_x, 2)
                i_yy += pow(j_y, 2)
                i_xy += (j_x * j_y)

                j_xh += (j_x * (w_halo[cell].h - w_c[i].h))
                j_yh += (j_y * (w_halo[cell].h - w_c[i].h))

                j_xhu += (j_x * (w_halo[cell].hu - w_c[i].hu))
                j_yhu += (j_y * (w_halo[cell].hu - w_c[i].hu))

                j_xhv += (j_x * (w_halo[cell].hv - w_c[i].hv))
                j_yhv += (j_y * (w_halo[cell].hv - w_c[i].hv))

                j_xhB1 += (j_x * (w_halo[cell].hB1 - w_c[i].hB1))
                j_yhB1 += (j_y * (w_halo[cell].hB1 - w_c[i].hB1))

                j_xhB2 += (j_x * (w_halo[cell].hB2 - w_c[i].hB2))
                j_yhB2 += (j_y * (w_halo[cell].hB2 - w_c[i].hB2))


                j_xz += (j_x * (w_halo[cell].Z - w_c[i].Z))
                j_yz += (j_y * (w_halo[cell].Z - w_c[i].Z))
        
        for k in range(3):
            nod = nodeidc[i][k]
                   
            if namen[nod] != 0 and namen[nod] != 10:
                for j in range(2):
                    center = centerg[nod][j]
                    cell = np.int(centerg[nod][j][2])
                    
                    j_x = center[0] - centerc[i][0]
                    j_y = center[1] - centerc[i][1]
                    i_xx += pow(j_x, 2)
                    i_yy += pow(j_y, 2)
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

#TODO
@njit(fastmath=True)
def face_gradient(wnd, Points, cellidf, nodeidf, centerg, namef, halofid, centerc, 
                  centerh, vertexn, w_xface, w_yface):
    
    nbface = len(cellidf)
    for i in range(nbface):

        c_left = cellidf[i][0]
        i_1 = nodeidf[i][0]
        i_2 = nodeidf[i][1]       
        
        xy_1 = vertexn[i_1][0:2]
        xy_2 = vertexn[i_2][0:2]
        
        v_1 = centerc[c_left]
        v_2 = np.array([Points[i][3]['x'], Points[i][3]['y']])

        f_1 = v_1 - xy_1
        f_2 = xy_2 - v_1
        f_3 = v_2 - xy_2
        f_4 = xy_1 - v_2
        
        
        hi1 = wnd[i]["h"][0]
        hi2 = wnd[i]["h"][1]
        hv1 = wnd[i]["h"][2]
        hv2 = wnd[i]["h"][3]
        
        ui1 = wnd[i]["hu"][0]/wnd[i]["h"][0]
        ui2 = wnd[i]["hu"][1]/wnd[i]["h"][1]
        uv1 = wnd[i]["hu"][2]/wnd[i]["h"][2]
        uv2 = wnd[i]["hu"][3]/wnd[i]["h"][3]
        
        vi1 = wnd[i]["hv"][0]/wnd[i]["h"][0]
        vi2 = wnd[i]["hv"][1]/wnd[i]["h"][1]
        vv1 = wnd[i]["hv"][2]/wnd[i]["h"][2]
        vv2 = wnd[i]["hv"][3]/wnd[i]["h"][3]
        
        zi1 = wnd[i]["Z"][0]
        zi2 = wnd[i]["Z"][1]
        zv1 = wnd[i]["Z"][2]
        zv2 = wnd[i]["Z"][3]
        
        a_i = 0.5 *((xy_2[0] - xy_1[0]) * (v_2[1]-v_1[1]) + (v_1[0]-v_2[0]) * (xy_2[1] - xy_1[1]))
   
        
        w_xface["hu"][i] = 1/(2*a_i)*((ui1 + uv1)*f_1[1] + (uv1 + ui2)*f_2[1] + (ui2 + uv2)*f_3[1] + (uv2 + ui1)*f_4[1])
        w_yface["hu"][i] = -1/(2*a_i)*((ui1 + uv1)*f_1[0] + (uv1 + ui2)*f_2[0] + (ui2 + uv2)*f_3[0] + (uv2 + ui1)*f_4[0])
             
        w_xface["hv"][i] = 1/(2*a_i)*((vi1 + vv1)*f_1[1] + (vv1 + vi2)*f_2[1] + (vi2 + vv2)*f_3[1] + (vv2 + vi1)*f_4[1])
        w_yface["hv"][i] = -1/(2*a_i)*((vi1 + vv1)*f_1[0] + (vv1 + vi2)*f_2[0] + (vi2 + vv2)*f_3[0] + (vv2 + vi1)*f_4[0])
        
        w_xface["h"][i] = 1/(2*a_i)*((hi1 + hv1)*f_1[1] + (hv1 + hi2)*f_2[1] + (hi2 + hv2)*f_3[1] + (hv2 + hi1)*f_4[1])
        w_yface["h"][i] = -1/(2*a_i)*((hi1 + hv1)*f_1[0] + (hv1 + hi2)*f_2[0] + (hi2 + hv2)*f_3[0] + (hv2 + hi1)*f_4[0])
        
        w_xface["Z"][i] = 1/(2*a_i)*((zi1 + zv1)*f_1[1] + (zv1 + zi2)*f_2[1] + (zi2 + zv2)*f_3[1] + (zv2 + zi1)*f_4[1])
        w_yface["Z"][i] = -1/(2*a_i)*((zi1 + zv1)*f_1[0] + (zv1 + zi2)*f_2[0] + (zi2 + zv2)*f_3[0] + (zv2 + zi1)*f_4[0])

    return w_xface, w_yface

@njit(fastmath=True)
def barthlimiter(w_c, w_ghost, w_halo, w_x, w_y, psi, cellid, faceid, namef, halofid, centerc, centerf):
    var = "h"
    nbelement = len(w_c)
    psi[:] = 2/3

    for i in range(nbelement):
        w_max = w_c[var][i]
        w_min = w_c[var][i]

        for j in range(3):
            face = faceid[i][j]
            if namef[face] == 0 and cellid[face][1] != -1:
                w_max = max(w_max, w_c[var][cellid[face][0]], w_c[var][cellid[face][1]])
                w_min = min(w_min, w_c[var][cellid[face][0]], w_c[var][cellid[face][1]])
            elif namef[face] != 10:
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
def update(w_c, wnew, dtime, rez, src, corio, vol):
    nbelement = len(w_c)

    for i in range(nbelement):

        wnew.h[i]   = w_c.h[i]   + dtime  *  ((rez["h"][i]   + src["h"][i])/vol[i]   + corio["h"][i])
        wnew.hu[i]  = w_c.hu[i]  + dtime  *  ((rez["hu"][i]  + src["hu"][i])/vol[i]  + corio["hu"][i])
        wnew.hv[i]  = w_c.hv[i]  + dtime  *  ((rez["hv"][i]  + src["hv"][i])/vol[i]  + corio["hv"][i])
        wnew.hB1[i] = w_c.hB1[i] + dtime  *  ((rez["hB1"][i] + src["hB1"][i])/vol[i] + corio["hB1"][i])
        wnew.hB2[i] = w_c.hB2[i] + dtime  *  ((rez["hB2"][i] + src["hB2"][i])/vol[i] + corio["hB2"][i])
        wnew.Z[i]   = w_c.Z[i]   + dtime  *  ((rez["Z"][i]   + src["Z"][i])/vol[i]   + corio["Z"][i])
        

    return wnew

@njit(fastmath=True)
def time_step(w_c, cfl, normal, mesure, volume, faceid, grav, term_convective):
    nbelement =  len(faceid)
    dt_c = 1e6
    kl=1

    for i in range(nbelement):
        velson = np.sqrt(grav*w_c[i].h)
        lam = 0
        if term_convective == "on" :
            for j in range(3):
                norm = normal[faceid[i][j]]
                mesn = mesure[faceid[i][j]]
                u_n = np.fabs(w_c.hu[i]/w_c.h[i]*norm[0] + w_c.hv[i]/w_c.h[i]*norm[1])
                u_n = u_n/mesn
                B_n = np.fabs(w_c.hB1[i]/w_c.h[i]*norm[0] + w_c.hB2[i]/w_c.h[i]*norm[1])
                B_n = B_n/mesn
                wb  = np.sqrt(B_n**2 + grav*w_c[i].h)
                lam1 = u_n - wb
                lam2 = u_n - B_n
                lam3 = u_n + B_n
                lam4 = u_n + wb
                lam_convect = max(lam1, lam2, lam3, lam4)
                lam += lam_convect * mesn
                
        if kl == 10 :
            for j in range(3):
                norm = normal[faceid[i][j]]
                mesn = mesure[faceid[i][j]]
                u_n = np.fabs(w_c.hu[i]/w_c.h[i]*norm[0] + w_c.hv[i]/w_c.h[i]*norm[1])
                lam_convect = u_n/mesn + velson
                lam += lam_convect * mesn

        if lam != 0:                 
            dt_c = min(dt_c, cfl * volume[i]/lam)

    dtime = np.asarray(dt_c)#np.min(dt_c))
    print(dt_c)

    return dtime

@njit(fastmath=True)
def exact_solution(wexact, center, time, grav, choix):
    nbelement = len(center)
    h_m = 2.534
    u_m = 4.03
    h_1 = 5
    h_2 = 1.
    
    h_m = 2.907
    u_m = 1.846
    h_1 = 4.
    h_2 = 2.   
    x_0 = 0.
    c_o = np.sqrt(grav*h_1)
    s_s = np.sqrt(0.5*grav*h_m*(1.+(h_m/h_2)))

    x_1 = x_0 - c_o*time
    x_2 = x_0 + 2.*c_o*time - 3.*np.sqrt(grav*h_m)*time
    x_3 = s_s*time + x_0
#    xl = 40.
#    cmax = 1.
#    sigma = 0.25
#    mu = 0.1
    
    if choix == 9:
        c1 = 0.04
        c2 = 0.02
        alpha = np.pi/6
        x0 = -20
        y0 = -10
        M = .5

    for i in range(nbelement):
        xcent = center[i][0]
        ycent = center[i][1]
        if choix == 0:
            

            if xcent < x_1:
                wexact[i].h = h_1
                wexact[i].hu = 0
        
            elif xcent < x_2 and xcent >= x_1:
                wexact[i].h = 1/(9*grav) * (2*np.sqrt(grav*h_1) - (xcent-x_0)/time)**2
                wexact[i].hu = (2/3 * (np.sqrt(grav*h_1) + (xcent-x_0)/time))*wexact[i].h
        
            elif xcent < x_3 and xcent >= x_2:
                wexact[i].h = h_m
                wexact[i].hu = u_m*wexact[i].h
        
            elif xcent >= x_3:
                wexact[i].h = h_2
                wexact[i].hu = 0
        
        elif choix == 1:
           
            wexact[i].Z = 10 + (40*xcent/14000) + 10*np.sin(np.pi*(4*xcent/14000 - 0.5))
            wexact[i].h = 64.5 - wexact[i].Z - 4*np.sin(np.pi*(4*time/86400 + 0.5))
            wexact[i].hu = (xcent - 14000)*np.pi/(5400*wexact[i].h)*np.cos(np.pi*(4*time/86400 + 0.5))
            
                
        elif choix == 3:

            if np.fabs(xcent - 1500/2) <= 1500/8 :
                wexact[i].Z = 8
            wexact[i].h = 20 - wexact[i].Z - 4*np.sin(np.pi*(4*time/86400 + 0.5))
            wexact[i].hu = (xcent - 1500)*np.pi/(5400*wexact[i].h) * np.cos(np.pi*(4*time/86400 + 0.5))
            
            wexact[i].hu = wexact[i].hu * wexact[i].h
            wexact[i].hv = wexact[i].hv * wexact[i].h
      
        elif choix == 9:

            f =  -c2*((xcent - x0 - M*time*np.cos(alpha))**2 + (ycent - y0 - M*time*np.sin(alpha))**2)
          
            wexact[i].h = 1 - ((c1**2)/(4*grav*c2))*np.exp(2*f)
            wexact[i].hu = M * np.cos(alpha) + c1*(ycent - y0 - M*time*np.sin(alpha))*np.exp(f)
            wexact[i].hv = M * np.sin(alpha) - c1*(xcent - x0 - M*time*np.cos(alpha))*np.exp(f)

            wexact[i].hu = wexact[i].hu * wexact[i].h
            wexact[i].hv = wexact[i].hv * wexact[i].h
        
#            xcent = center[i][0]
#            ycent = center[i][1]
#            hexact[i] = cmax/(1+(4*mu*time)/sigma**2)*np.exp(-(xcent**2+ycent**2)/(sigma**2 + 4*mu*time))


    return wexact

#from numba.typed import Dict

#TODO initialisation
@njit(fastmath=True)
def initialisation(w_c, center, choix):

    nbelements = len(center)

    if choix == 0: #Dam break
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            if (xcent <= 0.):
                w_c.h[i] = 4.
                w_c.hu[i] = 0.
                w_c.Z[i] = 0.
            else:
                w_c.h[i] = 2.
                w_c.Z[i] = 0.   
                w_c.hu[i] = 0
            
            w_c.hv[i] = 0.
            w_c.hB1[i] = 0.
            w_c.hB2[i] = 0.
    

    elif choix == 1: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            if xcent <= 0 :     
                w_c.h[i]  = 1.
                w_c.hu[i] = 0.
                w_c.hv[i] = 0.
                w_c.Z[i]  = 0.
                w_c.hB1[i] = 1.
                w_c.hB2[i] = 0.
            else:
                w_c.h[i]  = 2.
                w_c.Z[i]  = 0.   
                w_c.hu[i] = 0.
                w_c.hv[i] = 0.
                w_c.hB1[i] = 1.
                w_c.hB2[i] = 2.
            

            
                        
    elif choix == 2: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            if ycent <= 0 :     
                w_c.h[i]  = 1.
                w_c.hu[i] = 4.5
                w_c.hv[i] = 0.
                w_c.Z[i]  = 0.
                w_c.hB1[i] = 2.
                w_c.hB2[i] = 0.
            else:
                w_c.h[i]  = 2.
                w_c.Z[i]  = 0.   
                w_c.hu[i] = 11.
                w_c.hv[i] = 0.
                w_c.hB1[i] = 1.
                w_c.hB2[i] = 0.


    elif choix == 3: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            if np.sqrt(xcent**2 + ycent**2) <= 0.1 :     
                w_c.h[i]  = 10.
                w_c.hu[i] = -10.*ycent
                w_c.hv[i] = 10.*xcent
                w_c.Z[i]  = 0.
                w_c.hB1[i] = 1.
                w_c.hB2[i] = 0.
            else:
                w_c.h[i]  = 1.
                w_c.Z[i]  = 0.   
                w_c.hu[i] = 0.
                w_c.hv[i] = 0.
                w_c.hB1[i] = 1.
                w_c.hB2[i] = 0.


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
            w_c.hB1[i] = 0.0
            w_c.hB2[i] = 0.  



    return w_c



#TODO ghost
@njit
def ghost_value(w_c, w_ghost, cellid, name, normal, mesure, time, center, centerf,cellnid, w_x, w_y):


    nbface = len(cellid)

    for i in range(nbface):
        w_ghost[i] = w_c[cellid[i][0]]
        #xcent = center[i][0]
        #ycent = center[i][1]
        xc    = center[cellid[i][0]][0]
        yc    = center[cellid[i][0]][1]
        xe    = centerf[i][0]
        ye    = centerf[i][1]

        if (name[i] == 1 or name[i] == 2):         
          
            w_ghost[i].hu = w_c[cellid[i][0]].hu
            w_ghost[i].hv = w_c[cellid[i][0]].hv
            w_ghost[i].h  = w_c[cellid[i][0]].h 
            w_ghost[i].Z  = w_c[cellid[i][0]].Z
            w_ghost[i].hB1 = w_c[cellid[i][0]].hB1 
            w_ghost[i].hB2 = w_c[cellid[i][0]].hB2 
            
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
  
            w_ghost[i].h = w_c[cellid[i][0]].h
            w_ghost[i].Z = w_c[cellid[i][0]].Z
            w_ghost[i].hB1 = w_c[cellid[i][0]].hB1
            w_ghost[i].hB2 = w_c[cellid[i][0]].hB2
            w_ghost[i].hu = w_c[cellid[i][0]].h * u_g
            w_ghost[i].hv = w_c[cellid[i][0]].h * v_g
            
 


    return w_ghost



def save_paraview_results(w_c, wexact, niter, miter, time, dtime, rank, size, vol, cells, nodes, w_x, w_y):


    elements = {"triangle": cells}
    points = []
    for i in nodes:
        points.append([i[0], i[1], i[2]])
    
    cells  = np.array(cells)
    points = np.array(points)
    nbelement = len(w_c)
    Error_h = np.zeros(len(w_c))
    Error_u = np.zeros(len(w_c))
    he = np.zeros(len(w_c))
    ue = np.zeros(len(w_c))

#    for i in range(nbelement):
#
#        Error_h[i] = np.fabs(w_c["h"][i] - wexact["h"][i])*vol[i]
#        Error_u[i] = np.fabs(w_c["hu"][i]/w_c["h"][i] - wexact["hu"][i]/wexact["h"][i])*vol[i]
#        he[i]=np.fabs(wexact["h"][i])*vol[i]
#        ue[i]=np.fabs(wexact["hu"][i]/wexact["h"][i])*vol[i]     

#    ErrorL2_h =np.linalg.norm(Error_h,ord=1)/np.linalg.norm(he,ord=1)
#    ErrorL2_u =np.linalg.norm(Error_u,ord=1)/np.linalg.norm(ue,ord=1)
    data = {"h" : w_c["h"], "u" : w_c["hu"]/w_c["h"], "v": w_c["hv"]/w_c["h"], "B1": w_c["hB1"]/w_c["h"], "B2": w_c["hB2"]/w_c["h"], "Z": w_c["Z"], "h+Z": w_c["h"] + w_c["Z"], "uexact": wexact["hu"]/wexact["h"], "vexact": wexact["hv"]/wexact["h"],
            "hexact": wexact["h"], "divhB": w_x["hB1"] + w_y["hB2"]}
   
    if len(w_c) == len(cells):
        data = {"h": data, "u":data, "v": data, "B1": data, "B2": data, "Z":data, "h+Z":data, "uexact":data, "vexact":data,
                "hexact":data, "divhB":data}

    #maxhZ = max(w_c["h"]+w_c["Z"])
    maxhZ = max(w_x["hB1"] + w_y["hB2"])
    minhZ = min(np.fabs(w_x["hB1"] + w_y["hB2"]))
    maxu = max(w_c["hu"]/w_c["h"])
    
    integral_maxhZ = np.zeros(1)
    integral_minhZ = np.zeros(1)
    integral_maxu = np.zeros(1, dtype=np.double)

    COMM.Reduce(maxhZ, integral_maxhZ, MPI.MAX, 0)
    COMM.Reduce(minhZ, integral_minhZ, MPI.MIN, 0)
    COMM.Reduce(maxu, integral_maxu, MPI.MAX, 0)
    if rank == 0:
        print(" **************************** Computing ****************************")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Saving Results $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print("Iteration = ", niter, "time = ", np.float32(time), "time step = ", np.float32(dtime))
        #print("max h+Z =", np.float32(integral_maxhZ[0]), "max u =", np.float32(integral_maxu[0]))
        print("max divhB =", np.float32(integral_maxhZ[0]), "max u =", np.float32(integral_maxu[0]))
        print("min divhB =", np.float32(integral_minhZ[0]))
        #print("Error_L2 for h =", ErrorL2_h, "Error_L2 for u =", ErrorL2_u)
    
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
                text_file.write("<PDataArray type=\"Float64\" Name=\"Z\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"h+Z\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"uexact\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"vexact\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"hexact\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"divhB\" format=\"binary\"/>\n")
                text_file.write("</PCellData>\n")
            else:
                text_file.write("<PPointData Scalars=\"h\">\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"h\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"u\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"v\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"B1\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"B2\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"Z\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"h+Z\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"uexact\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"vexact\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"hexact\" format=\"binary\"/>\n")
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
