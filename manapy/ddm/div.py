#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:26:53 2020

@author: kissami
"""
from mpi4py import MPI
import numpy as np
from numba import njit
import manapy.ddm as ddm


COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

#TODO mhd_div
@njit(fastmath=True)
def compute_div(flux, fleft, fright,w_l, w_r, normal, mesure):

            norm = normal/mesure             
            Bnleft  = w_l["hB1"]*norm[0] + w_l["hB2"]*norm[1]	
            Bnright = w_r["hB1"]*norm[0] + w_r["hB2"]*norm[1]	    
            Bn = (Bnleft + Bnright)/2
            flux["h"]   = 0.0
            flux["hu"]  = 0.0
            flux["hv"]  = 0.0
            flux["hB1"] = 0.0
            flux["hB2"] = 0.0
            flux["hP"]   = 0.0
            flux["Ch"]   = Bn*mesure
            flux["Z"]   = 0
    

            return flux
#TODO
@njit(fastmath=True)
def explicitscheme_div(w_c, w_p, w_x, w_y, psi, w_ghost, w_halo, wx_halo, wy_halo, psi_halo, cellidf, faceidc,
                              nodeidc, centerc, cellidc, centerh, mesuref, centerf, normal, halofid,
                              name, ghostcenterf, cellidn, shift, mystruct, order, SCHEME,  grav):

    rezidus_div = np.zeros(len(w_c), dtype=mystruct)
    w_l  = np.zeros(1, dtype=mystruct)[0]
    w_r  = np.zeros(1, dtype=mystruct)[0]
    w_ln = np.zeros(1, dtype=mystruct)
    w_rn = np.zeros(1, dtype=mystruct)
    nbface = len(cellidf)
    
    flx    = np.zeros(1, dtype=mystruct)[0]
    fleft  = np.zeros(1, dtype=mystruct)[0]
    fright = np.zeros(1, dtype=mystruct)[0]

    if order == 1:

        for i in range(nbface):
            
            w_l = w_c[cellidf[i][0]]
            norm = normal[i]
            mesu = mesuref[i]           
            if name[i] == 0:
                w_r = w_c[cellidf[i][1]]
                flx = compute_div(flx, fleft, fright, w_l, w_r, norm, mesu)
                 
                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)
                rezidus_div[cellidf[i][1]] = ddm.add(rezidus_div[cellidf[i][1]], flx)

                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)
                rezidus_div[cellidf[i][1]] = ddm.add(rezidus_div[cellidf[i][1]], flx)

            elif name[i] == 5 or name[i] == 6 or name[i] == 7 or name[i] == 8:

                w_r = w_c[cellidf[i][1]]
                flx = compute_div(flx, fleft, fright, w_l, w_r, norm, mesu)
                                                  
                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)
                
            elif name[i] == 10:
                w_r = w_halo[halofid[i]]
                flx = compute_div(flx, fleft, fright, w_l, w_r, norm, mesu)

                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)

            else:
                w_r = w_ghost[i]
                flx = compute_div(flx, fleft, fright, w_l, w_r, norm, mesu)                
                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)


    elif order == 2:

        for i in range(nbface):
           
            norm = normal[i]
            mesu = mesuref[i]

            if name[i] == 0:
                w_l = w_c[cellidf[i][0]]
                w_r = w_c[cellidf[i][1]]
                
                center_left = centerc[cellidf[i][0]]
                center_right = centerc[cellidf[i][1]]

                w_x_left = w_x[cellidf[i][0]]
                w_y_left = w_y[cellidf[i][0]]
                psi_left = psi[cellidf[i][0]]

                w_x_right = w_x[cellidf[i][1]]
                w_y_right = w_y[cellidf[i][1]]
                psi_right = psi[cellidf[i][1]]

                r_l = np.array([centerf[i][0] - center_left[0], centerf[i][1] - center_left[1]])
                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])

                r_r = np.array([centerf[i][0] - center_right[0], centerf[i][1] - center_right[1]])

                w_rn["hB1"][0] = w_r["hB1"] + psi_right * (w_x_right["hB1"] * r_r[0] + w_y_right["hB1"] * r_r[1])
                w_rn["hB2"][0] = w_r["hB2"] + psi_right * (w_x_right["hB2"] * r_r[0] + w_y_right["hB2"] * r_r[1])
          
                flx = compute_div(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu)
                                                        
                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)
                rezidus_div[cellidf[i][1]] = ddm.add(rezidus_div[cellidf[i][1]], flx)

    
            elif name[i] == 5 or name[i] == 6 or name[i] == 7 or name[i] == 8:              

                w_l = w_c[cellidf[i][0]]
                w_r = w_c[cellidf[i][1]]
                
                center_left = centerc[cellidf[i][0]]
                center_right = centerc[cellidf[i][1]]

                w_x_left = w_x[cellidf[i][0]]
                w_y_left = w_y[cellidf[i][0]]
                psi_left = psi[cellidf[i][0]]

                w_x_right = w_x[cellidf[i][1]]
                w_y_right = w_y[cellidf[i][1]]
                psi_right = psi[cellidf[i][1]]

                r_l = np.array([centerf[i][0] - center_left[0], centerf[i][1] - center_left[1]])

                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])

                r_r = np.array([centerf[i][0] - center_right[0] - shift[cellidf[i][1]][0],
                                centerf[i][1] - center_right[1] - shift[cellidf[i][1]][1]])
                w_rn["hB1"][0] = w_r["hB1"] + psi_right * (w_x_right["hB1"] * r_r[0] + w_y_right["hB1"] * r_r[1])
                w_rn["hB2"][0] = w_r["hB2"] + psi_right * (w_x_right["hB2"] * r_r[0] + w_y_right["hB2"] * r_r[1])

                flx = compute_div(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu)
                             
                
                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)
                #rezidus_div[cellidf[i][1]] = ddm.add(rezidus_div[cellidf[i][1]], flx)

            elif name[i] == 10 and SIZE > 1:
                w_l = w_c[cellidf[i][0]]
                w_r = w_halo[halofid[i]]

                center_left = centerc[cellidf[i][0]]
                center_right = centerh[halofid[i]]

                w_x_left = w_x[cellidf[i][0]]
                w_y_left = w_y[cellidf[i][0]]
                psi_left = psi[cellidf[i][0]]

                w_x_right = wx_halo[halofid[i]]
                w_y_right = wy_halo[halofid[i]]
                psi_right = psi_halo[halofid[i]]

                r_l = np.array([centerf[i][0] - center_left[0], centerf[i][1] - center_left[1]])
                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])

                r_r = np.array([centerf[i][0] - centerh[halofid[i]][0],
                                centerf[i][1] - centerh[halofid[i]][1]])
 
                w_rn["hB1"][0] = w_r["hB1"] + psi_right * (w_x_right["hB1"] * r_r[0] + w_y_right["hB1"] * r_r[1])
                w_rn["hB2"][0] = w_r["hB2"] + psi_right * (w_x_right["hB2"] * r_r[0] + w_y_right["hB2"] * r_r[1])
 
                if SCHEME=="SRNH":
                    flx = compute_div(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu)

                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)

            else:
                w_l = w_c[cellidf[i][0]]
                w_r = w_ghost[i]


                center_left = centerc[cellidf[i][0]]
                center_right = ghostcenterf[i]

                w_x_left = w_x[cellidf[i][0]]
                w_y_left = w_y[cellidf[i][0]]
                psi_left = psi[cellidf[i][0]]

                r_l = np.array([centerf[i][0] - center_left[0], centerf[i][1] - center_left[1]])
                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])

                flx = compute_div(flx, fleft, fright, w_ln[0], w_r, norm, mesu)                                     
                rezidus_div[cellidf[i][0]] = ddm.minus(rezidus_div[cellidf[i][0]], flx)
                
    return rezidus_div


def calculate_div(w_c, rez, vol):
    nbelement = len(w_c)
    div = np.zeros(len(w_c))

    for i in range(nbelement):


        div[i] = rez["Ch"][i]/vol[i]         
        

    return div





@njit(fastmath=True)
def source_GLM_correction(w_c, dtime, mystruct):
    coriolis = np.zeros(len(w_c), dtype=mystruct)
    nbelement = len(w_c)
    
    for i in range(nbelement):

        coriolis["h"][i] = 0.
        coriolis["hu"][i] = 0#f_c*w_c["hv"][i]
        coriolis["hv"][i] = 0#-f_c*w_c["hu"][i]
        coriolis["hB1"][i] = 0.
        coriolis["hB2"][i] = 0.
        coriolis["Z"][i] = 0.
    
    return coriolis    
