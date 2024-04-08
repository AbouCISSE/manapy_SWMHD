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
from numpy import zeros, array


COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

@njit
def compute_flux_fvc(flux, w_p, normal, grav, mesure):

    hu = w_p["hu"]
    hv = w_p["hv"]
    hB1 = w_p["hB1"]
    hB2 = w_p["hB2"]
    h  = w_p["h"]
   
    flux["h"]  = hu*normal[0] + hv*normal[1]
    flux["hu"] = (hu*hu/h + 0.5*grav*h*h - hB1*hB1/h)*normal[0] + (hu*hv/h - hB1*hB2/h)*normal[1]
    flux["hv"] = (hu*hv/h - hB1*hB2/h)*normal[0] + (hv*hv/h + 0.5*grav*h*h - hB2*hB2/h )*normal[1]
    flux["hB1"] = (hv*hB1/h - hu*hB2/h)*normal[1]
    flux["hB2"] = (hu*hB2/h - hv*hB1/h)*normal[0]
    flux["Z"]  = 0
    
    return flux
 
#TODO mhd_10
@njit(fastmath=True)
def compute_flux_shallow_srnh_mhd(flux, fleft, fright, w_l, w_r, normal, mesure, grav, w_pl, w_pr):


    ninv =np.zeros(2)
    norm = normal/mesure
    ninv[0] = -1*normal[1]
    ninv[1] = normal[0]
    w_dif =np.zeros(6)
    
    if (w_l["h"] < 1e-9 or w_r["h"] < 1e-9) : 
        print("hl is zero ", w_l["h"])
 


    u_h = (w_l["hu"] / w_l["h"] * np.sqrt(w_l["h"])
       + w_r["hu"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
    v_h = (w_l["hv"] / w_l["h"] * np.sqrt(w_l["h"])
           + w_r["hv"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
           
    B1_h = (w_l["hB1"] / w_l["h"] * np.sqrt(w_l["h"])
       + w_r["hB1"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))

    B2_h = (w_l["hB2"] / w_l["h"] * np.sqrt(w_l["h"])
           + w_r["hB2"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
    
 
    un_h = u_h*normal[0] + v_h*normal[1]
    un_h = un_h / mesure
    vn_h = u_h*ninv[0] + v_h*ninv[1]
    vn_h = vn_h / mesure
    
    B1n_h = B1_h*normal[0] + B2_h*normal[1]
    B1n_h = B1n_h / mesure
    B2n_h = B1_h*ninv[0] + B2_h*ninv[1]
    B2n_h = B2n_h / mesure

    hroe  = (w_l["h"]+w_r["h"])/2
    uroe  = un_h
    vroe  = vn_h
    B1roe = B1n_h
    B2roe = B2n_h

    
    
    
    uleft = w_l["hu"]*normal[0] + w_l["hv"]*normal[1]
    uleft = uleft / mesure
    
    vleft = w_l["hu"]*ninv[0] + w_l["hv"]*ninv[1]
    vleft = vleft / mesure

    uright = w_r["hu"]*normal[0] + w_r["hv"]*normal[1]
    uright = uright / mesure
    
    vright = w_r["hu"]*ninv[0] + w_r["hv"]*ninv[1]
    vright = vright / mesure

    B1left = w_l["hB1"]*normal[0] + w_l["hB2"]*normal[1]
    B1left = B1left / mesure
    
    B2left = w_l["hB1"]*ninv[0] + w_l["hB2"]*ninv[1]
    B2left = B2left / mesure

    B1right = w_r["hB1"]*normal[0] + w_r["hB2"]*normal[1]
    B1right = B1right / mesure
    
    B2right = w_r["hB1"]*ninv[0] + w_r["hB2"]*ninv[1]
    B2right = B2right / mesure
    
    B1i = (w_pl["hB1"]*normal[0] + w_pl["hB2"]*normal[1]) / mesure
    B1j = (w_pr["hB1"]*normal[0] + w_pr["hB2"]*normal[1])/ mesure
    
    hs = np.sqrt(w_l["h"]) + np.sqrt(w_r["h"])
    hd = (w_r["h"] + w_l["h"])/2
    B1h = w_l["hB1"] / w_l["h"] * np.sqrt(w_l["h"])+ w_r["hB1"] / w_r["h"] * np.sqrt(w_r["h"])
    B2h = w_l["hB2"] / w_l["h"] * np.sqrt(w_l["h"])+ w_r["hB2"] / w_r["h"] * np.sqrt(w_r["h"]) 
    B1i = (w_pl["hB1"]*normal[0] + w_pl["hB2"]*normal[1]) / mesure
    B1j = (w_pr["hB1"]*normal[0] + w_pr["hB2"]*normal[1])/ mesure     
                       
    absy = (B1i + B1j)/2  

    
    w_lrh = (w_l["h"]  + w_r["h"])/2
    w_lrhu = (uleft + uright)/2
    w_lrhv = (vleft + vright)/2
    w_lrhB1 = (B1left + B1right)/2
    w_lrhB2 = (B2left + B2right)/2
    w_lrz = (w_l["Z"] + w_r["Z"])/2
    w_lrhP = (w_l["hP"]  + w_r["hP"])/2
    
                        
    w_dif[0] = w_r["h"] - w_l["h"]
    w_dif[1] = uright   - uleft
    w_dif[2] = vright   - vleft
    w_dif[3] = B1right  - B1left
    w_dif[4] = B2right  - B2left
    w_dif[5] = w_r["Z"] - w_l["Z"]
    
    signA= np.zeros((6, 6))
        
    sound = np.sqrt(grav * hroe)
    

    
    w=np.sqrt(B1roe * B1roe + grav * hroe)

    
    lambda1 = uroe - w
    #print("uroe =",uroe)
    #print("w =",w)
    
    
    lambda2 = uroe - B1roe
    
    lambda3 = uroe + B1roe
    
    lambda4 = uroe + w
    
    cpsi = max(np.fabs(lambda1) , np.fabs(lambda4))
    

 
    eps=1e-15
    if np.fabs(lambda1) < eps :
            s1  = 0.
            pi1 = 0.
            #print("lambda1=", lambda1)
    else :
            s1  = lambda1/np.fabs(lambda1)
            pi1 = s1/lambda1

            
    if np.fabs(lambda2) < eps :
            s2  = 0.
            pi2 = 0.
            #print("lambda2=", lambda2)
            #print("Champroe=",B1roe-uroe)
    else :
            s2  = lambda2/np.fabs(lambda2)  
            pi2 = 1/np.fabs(lambda2)        
            
    if np.fabs(lambda3) < eps :
            s3  = 0.
            pi3 = 0.
           # print("lambda3=", lambda3)
    else :
            s3  = lambda3/np.fabs(lambda3)
            pi3 = 1/np.fabs(lambda3)
            
    if np.fabs(lambda4) < eps :
            s4  = 0.
            pi4 = 0.
            #print("lambda4=", lambda4)
    else :
            s4  = lambda4/np.fabs(lambda4)
            pi4 = 1/np.fabs(lambda4)
            

    gamma1 = vroe + B2roe
    gamma2 = vroe - B2roe
         
    sigma1 = vroe*(s1*lambda4  - s4*lambda1)  - w*(s2*gamma1 + s3*gamma2)
    sigma2 = B2roe*(s1*lambda4 - s4*lambda1)  - w*(s2*gamma1 - s3*gamma2)
    
    if (np.fabs(lambda2) < eps and np.fabs(lambda3)) < eps  :

            mu1 = B1roe*vroe*pi1/w  - B1roe*vroe*pi4/w
            mu2 = B1roe*B2roe*pi1/w - B1roe*B2roe*pi4/w
            ann = 0 
    else :   
            mu1 = B1roe*vroe*pi1/w  - B1roe*vroe*pi4/w  - 0.5*(gamma1*pi2 - gamma2*pi3)  
            mu2 = B1roe*B2roe*pi1/w - B1roe*B2roe*pi4/w - 0.5*(gamma1*pi2 +  gamma2*pi3) 
            ann = 0.0

            
        
    #1ère colonne de la matrice A
    signA[0][0] = (s1*lambda4 - s4*lambda1)/(2*w)
    signA[1][0] = lambda1*lambda4*(s1-s4)/(2*w)
    signA[2][0] = sigma1/(2*w)
    signA[3][0] = 0.0
    signA[4][0] = sigma2/(2*w)
    signA[5][0] = 0.0
                
    #2ème colonne de la matrice A
    signA[0][1] = (s4 - s1)/(2*w)
    signA[1][1] = (s4*lambda4 - s1*lambda1)/(2*w)
    signA[2][1] = vroe*(s4 - s1)/(2*w)
    signA[3][1] = 0.0
    signA[4][1] = B2roe*(s4 - s1)/(2*w)
           
    signA[5][1] = 0.0
                
    #3ème colonne de la matrice A
    signA[0][2] = 0.0
    signA[1][2] = 0.0
    signA[2][2] = (s2 + s3)/2
    signA[3][2] = 0.0
    signA[4][2] = (s2 - s3)/2
    signA[5][2] = 0.0
                
    #4ème colonne de la matrice A
    signA[0][3] = ann*B1roe*(pi1 - pi4)/w
    signA[1][3] = ann*B1roe*(s1-s4)/w
    signA[2][3] = ann*mu1
    signA[3][3] = 0.0
    signA[4][3] = ann*mu2
    signA[5][3] = 0.0
                
    #5ème colonne de la matrice A
    signA[0][4] = 0.0
    signA[1][4] = 0.0
    signA[2][4] = (s2 - s3)/2
    signA[3][4] = 0.0
    signA[4][4] = (s2 + s3)/2
    signA[5][4] = 0.0
                
    #6ème colonne de la matrice A
    signA[0][5] = (sound**2)*(pi4 - pi1)/(2*w)
    signA[1][5] = (sound**2)*(s4-s1)/(2*w)
    signA[2][5] = (sound**2)*vroe*(pi4 - pi1)/(2*w)
    signA[3][5] = 0.0
    signA[4][5] = (sound**2)*B2roe*(pi4 - pi1)/(2*w)
    signA[5][5] = 0.0

    
        

    smmat=signA
    
    hnew  = 0.
    unew  = 0.
    vnew  = 0.
    B1new = 0.
    B2new = 0.
    znew  = 0.

    for i in range(6):
        hnew  += smmat[0][i] * w_dif[i]
        unew  += smmat[1][i] * w_dif[i]
        vnew  += smmat[2][i] * w_dif[i]
        B1new += smmat[3][i] * w_dif[i]
        B2new += smmat[4][i] * w_dif[i]
        znew  += smmat[5][i] * w_dif[i]



    Pnew  = cpsi * (B1right - B1left)
    u_h   = hnew/2
    u_hu  = unew/2
    u_hv  = vnew/2
    u_hP  = Pnew/2
    u_hB1 = B1new/2
    u_hB2 = B2new/2
    u_z   = znew/2

    w_lrh   = w_lrh   - u_h
    w_lrhu  = w_lrhu  - u_hu
    w_lrhv  = w_lrhv  - u_hv
    w_lrhP  = w_lrhP  - u_hP
    w_lrhB1 = w_lrhB1 - u_hB1
    w_lrhB2 = w_lrhB2 - u_hB2
    w_lrz   = w_lrz   - u_z
    
    
 
    
    w_lrhB1 = absy
    unew  = 0.
    vnew  = 0.
    B1new = 0.
    B2new = 0.
    
    
    unew = w_lrhu * normal[0]    - w_lrhv  * normal[1]
    unew = unew / mesure
    vnew = w_lrhu * normal[1]    + w_lrhv  * normal[0] 
    vnew = vnew / mesure

    B1new = w_lrhB1 * normal[0]  - w_lrhB2 * normal[1]
    B1new = B1new / mesure
    B2new = w_lrhB1 * normal[1]  + w_lrhB2 * normal[0] 
    B2new = B2new / mesure


    w_lrhu = unew
    w_lrhv = vnew
    
    
    w_lrhB1 = B1new
    w_lrhB2 = B2new

    q_s = normal[0] * unew  + normal[1] * vnew
    p_s = normal[0] * B1new + normal[1] * B2new
    

    flux["h"]   = q_s
    flux["hu"]  = q_s * w_lrhu/w_lrh + 0.5 * grav * w_lrh * w_lrh * normal[0] - p_s*w_lrhB1/w_lrh
    flux["hv"]  = q_s * w_lrhv/w_lrh + 0.5 * grav * w_lrh * w_lrh * normal[1] - p_s*w_lrhB2/w_lrh
    flux["hB1"] = (w_lrhv*w_lrhB1/w_lrh - w_lrhu*w_lrhB2/w_lrh ) * normal[1] 
    flux["hB2"] = (w_lrhu*w_lrhB2/w_lrh - w_lrhv*w_lrhB1/w_lrh ) * normal[0] 
    flux["hP"]   = 0.0
    flux["Ch"]   = 0.0
    flux["Z"]   = 0 





    return flux
 
 
 
#TODO mhd_10
@njit(fastmath=True)
def compute_flux_shallow_srnh_mhd_GLM(flux, fleft, fright, w_l, w_r, normal, mesure, grav, w_pl, w_pr, cpsi):


    ninv =np.zeros(2)
    norm = normal/mesure
    ninv[0] = -1*normal[1]
    ninv[1] = normal[0]
    w_dif =np.zeros(6)


    u_h = (w_l["hu"] / w_l["h"] * np.sqrt(w_l["h"])
       + w_r["hu"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
    v_h = (w_l["hv"] / w_l["h"] * np.sqrt(w_l["h"])
           + w_r["hv"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
           
    B1_h = (w_l["hB1"] / w_l["h"] * np.sqrt(w_l["h"])
       + w_r["hB1"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))

    B2_h = (w_l["hB2"] / w_l["h"] * np.sqrt(w_l["h"])
           + w_r["hB2"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
    
 
    un_h = u_h*normal[0] + v_h*normal[1]
    un_h = un_h / mesure
    vn_h = u_h*ninv[0] + v_h*ninv[1]
    vn_h = vn_h / mesure
    
    B1n_h = B1_h*normal[0] + B2_h*normal[1]
    B1n_h = B1n_h / mesure
    B2n_h = B1_h*ninv[0] + B2_h*ninv[1]
    B2n_h = B2n_h / mesure

    hroe  = (w_l["h"]+w_r["h"])/2
    uroe  = un_h
    vroe  = vn_h
    B1roe = B1n_h
    B2roe = B2n_h

    
    
    
    uleft = w_l["hu"]*normal[0] + w_l["hv"]*normal[1]
    uleft = uleft / mesure
    
    vleft = w_l["hu"]*ninv[0] + w_l["hv"]*ninv[1]
    vleft = vleft / mesure

    uright = w_r["hu"]*normal[0] + w_r["hv"]*normal[1]
    uright = uright / mesure
    
    vright = w_r["hu"]*ninv[0] + w_r["hv"]*ninv[1]
    vright = vright / mesure

    B1left = w_l["hB1"]*normal[0] + w_l["hB2"]*normal[1]
    B1left = B1left / mesure
    
    B2left = w_l["hB1"]*ninv[0] + w_l["hB2"]*ninv[1]
    B2left = B2left / mesure

    B1right = w_r["hB1"]*normal[0] + w_r["hB2"]*normal[1]
    B1right = B1right / mesure
    
    B2right = w_r["hB1"]*ninv[0] + w_r["hB2"]*ninv[1]
    B2right = B2right / mesure
    
    B1i = (w_pl["hB1"]*normal[0] + w_pl["hB2"]*normal[1]) / mesure
    B1j = (w_pr["hB1"]*normal[0] + w_pr["hB2"]*normal[1])/ mesure     
                       
    absy = (B1i + B1j)/2 

    
    w_lrh = (w_l["h"]  + w_r["h"])/2
    w_lrhu = (uleft + uright)/2
    w_lrhv = (vleft + vright)/2
    w_lrhB1 = (B1left + B1right)/2
    w_lrhB2 = (B2left + B2right)/2
    w_lrz = (w_l["Z"] + w_r["Z"])/2
    w_lrhP = (w_l["hP"]  + w_r["hP"])/2
    
                        
    w_dif[0] = w_r["h"] - w_l["h"]
    w_dif[1] = uright   - uleft
    w_dif[2] = vright   - vleft
    w_dif[3] = B1right  - B1left
    w_dif[4] = B2right  - B2left
    w_dif[5] = w_r["Z"] - w_l["Z"]
    
    signA= np.zeros((6, 6))
        
    sound = np.sqrt(grav * hroe)
    
    B1roe = absy/hroe
    
    w=np.sqrt(B1roe * B1roe + grav * hroe)

    
    lambda1 = uroe - w
    
    lambda2 = uroe - B1roe
    
    lambda3 = uroe + B1roe
    
    lambda4 = uroe + w
    
    #cpsi = max(np.fabs(lambda1) , np.fabs(lambda4))
        
    eps=1e-15
    if np.fabs(lambda1) < eps :
            s1  = 0.
            pi1 = 0.
    else :
            s1  = lambda1/np.fabs(lambda1)
            pi1 = s1/lambda1

            
    if np.fabs(lambda2) < eps :
            s2  = 0.
            pi2 = 0.

    else :
            s2  = lambda2/np.fabs(lambda2)  
            pi2 = 1/np.fabs(lambda2)        
            
    if np.fabs(lambda3) < eps :
            s3  = 0.
            pi3 = 0.
    else :
            s3  = lambda3/np.fabs(lambda3)
            pi3 = 1/np.fabs(lambda3)
            
    if np.fabs(lambda4) < eps :
            s4  = 0.
            pi4 = 0.
    else :
            s4  = lambda4/np.fabs(lambda4)
            pi4 = 1/np.fabs(lambda4)
            

    gamma1 = vroe + B2roe
    gamma2 = vroe - B2roe
         
    sigma1 = vroe*(s1*lambda4  - s4*lambda1)  - w*(s2*gamma1 + s3*gamma2)
    sigma2 = B2roe*(s1*lambda4 - s4*lambda1)  - w*(s2*gamma1 - s3*gamma2)
    
    if (np.fabs(lambda2) < eps and np.fabs(lambda3)) < eps  :

            mu1 = B1roe*vroe*pi1/w  - B1roe*vroe*pi4/w
            mu2 = B1roe*B2roe*pi1/w - B1roe*B2roe*pi4/w
            ann = 0          

    else :   
            mu1 = B1roe*vroe*pi1/w  - B1roe*vroe*pi4/w  - 0.5*(gamma1*pi2 - gamma2*pi3)  
            mu2 = B1roe*B2roe*pi1/w - B1roe*B2roe*pi4/w - 0.5*(gamma1*pi2 +  gamma2*pi3) 
            ann = 0.0

        
    #1ère colonne de la matrice A
    signA[0][0] = (s1*lambda4 - s4*lambda1)/(2*w)
    signA[1][0] = lambda1*lambda4*(s1-s4)/(2*w)
    signA[2][0] = sigma1/(2*w)
    signA[3][0] = 0.0
    signA[4][0] = sigma2/(2*w)
    signA[5][0] = 0.0
                
    #2ème colonne de la matrice A
    signA[0][1] = (s4 - s1)/(2*w)
    signA[1][1] = (s4*lambda4 - s1*lambda1)/(2*w)
    signA[2][1] = vroe*(s4 - s1)/(2*w)
    signA[3][1] = 0.0
    signA[4][1] = B2roe*(s4 - s1)/(2*w)
           
    signA[5][1] = 0.0
                
    #3ème colonne de la matrice A
    signA[0][2] = 0.0
    signA[1][2] = 0.0
    signA[2][2] = (s2 + s3)/2
    signA[3][2] = 0.0
    signA[4][2] = (s2 - s3)/2
    signA[5][2] = 0.0
                
    #4ème colonne de la matrice A
    signA[0][3] = ann*B1roe*(pi1 - pi4)/w
    signA[1][3] = ann*B1roe*(s1-s4)/w
    signA[2][3] = ann*mu1
    signA[3][3] = 0.0
    signA[4][3] = ann*mu2
    signA[5][3] = 0.0
                
    #5ème colonne de la matrice A
    signA[0][4] = 0.0
    signA[1][4] = 0.0
    signA[2][4] = (s2 - s3)/2
    signA[3][4] = 0.0
    signA[4][4] = (s2 + s3)/2
    signA[5][4] = 0.0
                
    #6ème colonne de la matrice A
    signA[0][5] = (sound**2)*(pi4 - pi1)/(2*w)
    signA[1][5] = (sound**2)*(s4-s1)/(2*w)
    signA[2][5] = (sound**2)*vroe*(pi4 - pi1)/(2*w)
    signA[3][5] = 0.0
    signA[4][5] = (sound**2)*B2roe*(pi4 - pi1)/(2*w)
    signA[5][5] = 0.0
        

    smmat=signA
    
    hnew  = 0.
    unew  = 0.
    vnew  = 0.
    B1new = 0.
    B2new = 0.
    znew  = 0.

    for i in range(6):
        hnew  += smmat[0][i] * w_dif[i]
        unew  += smmat[1][i] * w_dif[i]
        vnew  += smmat[2][i] * w_dif[i]
        B1new += smmat[3][i] * w_dif[i]
        B2new += smmat[4][i] * w_dif[i]
        znew  += smmat[5][i] * w_dif[i]


    Pnew  = cpsi * (B1right - B1left)
    u_h   = hnew/2
    u_hu  = unew/2
    u_hv  = vnew/2
    u_hP  = Pnew/2
    u_hB1 = B1new/2
    u_hB2 = B2new/2
    u_z   = znew/2

    w_lrh   = w_lrh   - u_h
    w_lrhu  = w_lrhu  - u_hu
    w_lrhv  = w_lrhv  - u_hv
    w_lrhP  = w_lrhP  - u_hP
    w_lrhB1 = absy#w_lrhB1 - u_hB1
    w_lrhB2 = w_lrhB2 - u_hB2
    w_lrz   = w_lrz   - u_z

    w_hP   = w_r["hP"] - w_l["hP"]
    mw_hB1 = w_lrhB1/2 - (1/(2*cpsi))*w_hP
    mhP  = (w_r["hP"] + w_l["hP"])/2 - cpsi*w_dif[3]/2
    

    unew  = 0.
    vnew  = 0.
    B1new = 0.
    B2new = 0.
    
    
    unew = w_lrhu * normal[0]    - w_lrhv  * normal[1]
    unew = unew / mesure
    vnew = w_lrhu * normal[1]    + w_lrhv  * normal[0] 
    vnew = vnew / mesure

    B1new = w_lrhB1 * normal[0]  - w_lrhB2 * normal[1]
    B1new = B1new / mesure
    B2new = w_lrhB1 * normal[1]  + w_lrhB2 * normal[0] 
    B2new = B2new / mesure


    w_lrhu = unew
    w_lrhv = vnew
    
    
    w_lrhB1 = B1new
    w_lrhB2 = B2new

    q_s = normal[0] * unew  + normal[1] * vnew
    p_s = normal[0] * B1new + normal[1] * B2new
    
    Flux_B1psi = mhP * norm[0]*mesure
    Flux_B2psi = mhP * norm[1]*mesure
    Flux_hPpsi = cpsi*cpsi*mw_hB1*mesure
                       
                       
            

    flux["h"]   = q_s
    flux["hu"]  = q_s * w_lrhu/w_lrh + 0.5 * grav * w_lrh * w_lrh * normal[0] - p_s*w_lrhB1/w_lrh
    flux["hv"]  = q_s * w_lrhv/w_lrh + 0.5 * grav * w_lrh * w_lrh * normal[1] - p_s*w_lrhB2/w_lrh
    flux["hB1"] = (w_lrhv*w_lrhB1/w_lrh - w_lrhu*w_lrhB2/w_lrh ) * normal[1]  + Flux_B1psi 
    flux["hB2"] = (w_lrhu*w_lrhB2/w_lrh - w_lrhv*w_lrhB1/w_lrh ) * normal[0]  + Flux_B2psi  
    flux["hP"]   = Flux_hPpsi
    flux["Ch"]   = cpsi
    flux["Z"]   = 0 





    return flux 
 




 
 
 
 
 
 
 



#TODO mhd_roe
@njit(fastmath=True)
def compute_flux_swmhd_LF(flux, fleft, fright, w_l, w_r, normal, mesure, grav):
        #print("Im Roe scheme")
            norm = normal/mesure
            ninv =np.zeros(2)
            ninv[0] = -1*norm[1]
            ninv[1] = norm[0]

                
            hl  = w_l["h"]    
            ul  = w_l["hu"] / w_l["h"] 
            vl  = w_l["hv"] / w_l["h"] 
            B1l = w_l["hB1"] / w_l["h"] 
            B2l = w_l["hB2"] / w_l["h"] 
            hr  =  w_r['h']   
            ur  = w_r["hu"] / w_r["h"] 
            vr  = w_r["hv"] / w_r["h"] 
            B1r = w_r["hB1"] / w_r["h"] 
            B2r = w_r["hB2"] / w_r["h"] 
            
            u_h = (w_l["hu"] / w_l["h"] * np.sqrt(w_l["h"])
	          + w_r["hu"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
            v_h = (w_l["hv"] / w_l["h"] * np.sqrt(w_l["h"])
	    	   + w_r["hv"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
            B1_h = (w_l["hB1"] / w_l["h"] * np.sqrt(w_l["h"])
	       + w_r["hB1"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
	    
            B2_h = (w_l["hB2"] / w_l["h"] * np.sqrt(w_l["h"])
		   + w_r["hB2"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))  
            un_h = u_h*normal[0] + v_h*normal[1] 
            un_h = un_h / mesure
            vn_h = u_h*ninv[0] + v_h*ninv[1]
            vn_h = vn_h / mesure
            B1n_h = B1_h*normal[0] + B2_h*normal[1]
            B1n_h = B1n_h / mesure
            B2n_h = B1_h*ninv[0] + B2_h*ninv[1]
            B2n_h = B2n_h / mesure
            hroe  = (w_l["h"]+w_r["h"])/2
            uroe  = un_h
            vroe  = vn_h
            B1roe = B1n_h
            B2roe = B2n_h
            n1     = norm[0]
            n2     = norm[1]
            
            Unl     = ul*n1  + vl*n2   
            Bnl     = B1l*n1 + B2l*n2
            wl      = np.sqrt(grav * hl + Bnl**2)
            Unr     = ur*n1  + vr*n2
            Bnr     = B1r*n1 + B2r*n2
            wr      = np.sqrt(grav * hr + Bnr**2)
                        
#            Utr     = vr*n1  - ur*n2 
#            Btr     = B2r*n1 - B1r*n2   
#            Utl     = vl*n1  - ul*n2 
#            Btl     = B2l*n1 - B1l*n2
                     
    

            lambda1l = np.fabs(Unl - wl)
            lambda2l = np.fabs(Unl - Bnl)
            lambda3l = np.fabs(Unl + Bnl)
            lambda4l = np.fabs(Unl + wl) 
            
            lambda1r = np.fabs(Unr - wr)
            lambda2r = np.fabs(Unr - Bnr)
            lambda3r = np.fabs(Unr + Bnr)
            lambda4r = np.fabs(Unr + wr)
            
            ll = max(lambda1l, lambda2l, lambda3l, lambda4l)
            lr = max(lambda1r, lambda2r, lambda3r, lambda4r)
            
            lambda_star = max(ll,lr)
    
            q_l = w_l["hu"] * norm[0] + w_l["hv"] * norm[1]
            m_l = w_l["hB1"] * norm[0] + w_l["hB2"] * norm[1]
            p_l = 0.5 * grav * w_l["h"]* w_l["h"]
            
            q_r = w_r["hu"] * norm[0] + w_r["hv"] * norm[1]
            m_r = w_r["hB1"] * norm[0] + w_r["hB2"] * norm[1]
            p_r = 0.5 * grav * w_r["h"] * w_r["h"]

            fleft["h"]  = q_l
            fleft["hu"]  = q_l * w_l["hu"] /w_l["h"]  + p_l*norm[0] - m_l * w_l['hB1'] /w_l["h"]
            fleft["hv"]  = q_l * w_l["hv"] /w_l["h"]  + p_l*norm[1] - m_l * w_l['hB2'] /w_l["h"]
            fleft["hB1"] = (w_l["hv"] * w_l['hB1'] /w_l["h"] - w_l["hu"] * w_l['hB2'] /w_l["h"] )*norm[1]
            fleft["hB2"] = (w_l["hu"] * w_l['hB2'] /w_l["h"] - w_l["hv"] * w_l['hB1'] /w_l["h"])* norm[0]
  
            fright["h"] = q_r
            fright["hu"]  = q_r * w_r["hu"] /w_r["h"]  + p_r*norm[0] - m_r * w_r['hB1'] /w_r["h"]
            fright["hv"]  = q_r * w_r["hv"] /w_r["h"]  + p_r*norm[1] - m_r * w_r['hB2'] /w_r["h"]
            fright["hB1"] = (w_r["hv"] * w_r['hB1'] /w_r["h"] - w_r["hu"] * w_r['hB2'] /w_r["h"] )*norm[1]
            fright["hB2"] = (w_r["hu"] * w_r['hB2'] /w_r["h"] - w_r["hv"] * w_r['hB1'] /w_r["h"])* norm[0]
        
            f_h = 0.5 * (fleft["h"] + fright["h"])       - 0.5*lambda_star*(w_r["h"]   - w_l["h"])
            f_hu = 0.5 * (fleft["hu"] + fright["hu"])    - 0.5*lambda_star*(w_r["hu"]  - w_l["hu"])
            f_hv = 0.5 * (fleft["hv"] + fright["hv"])    - 0.5*lambda_star*(w_r["hv"]  - w_l["hv"])
            f_hB1 = 0.5 * (fleft["hB1"] + fright["hB1"]) - 0.5*lambda_star*(w_r["hB1"] - w_l["hB1"])
            f_hB2 = 0.5 * (fleft["hB2"] + fright["hB2"]) - 0.5*lambda_star*( w_r["hB2"] - w_l["hB2"])
            
          
            
            flux["h"]   = f_h * mesure
            flux["hu"]  = f_hu * mesure
            flux["hv"]  = f_hv * mesure
            flux["hB1"] = f_hB1 * mesure
            flux["hB2"] = f_hB2 * mesure
            flux["hP"]  = 0
            flux['Z']   = 0

            return flux


#TODO mhd_roe
@njit(fastmath=True)
def compute_flux_swmhd_LF_GLM_correction(flux, fleft, fright, w_l, w_r, normal, mesure, grav, c_h, delta):

            norm = normal/mesure
            ninv =np.zeros(2)
            ninv[0] = -1*norm[1]
            ninv[1] = norm[0]            
            tol = 1e-9
            
                
            hl  = w_l["h"]    
            ul  = w_l["hu"] / w_l["h"] 
            vl  = w_l["hv"] / w_l["h"] 
            B1l = w_l["hB1"] / w_l["h"] 
            B2l = w_l["hB2"] / w_l["h"] 
            hr  =  w_r['h']   
            ur  = w_r["hu"] / w_r["h"] 
            vr  = w_r["hv"] / w_r["h"] 
            B1r = w_r["hB1"] / w_r["h"] 
            B2r = w_r["hB2"] / w_r["h"] 

            u_h = (w_l["hu"] / w_l["h"] * np.sqrt(w_l["h"])
	          + w_r["hu"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
            v_h = (w_l["hv"] / w_l["h"] * np.sqrt(w_l["h"])
	    	   + w_r["hv"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
            B1_h = (w_l["hB1"] / w_l["h"] * np.sqrt(w_l["h"])
	       + w_r["hB1"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))
	    
            B2_h = (w_l["hB2"] / w_l["h"] * np.sqrt(w_l["h"])
		   + w_r["hB2"] / w_r["h"] * np.sqrt(w_r["h"])) /(np.sqrt(w_l["h"]) + np.sqrt(w_r["h"]))  
            un_h = u_h*normal[0] + v_h*normal[1] 
            un_h = un_h / mesure
            vn_h = u_h*ninv[0] + v_h*ninv[1]
            vn_h = vn_h / mesure
            B1n_h = B1_h*normal[0] + B2_h*normal[1]
            B1n_h = B1n_h / mesure
            B2n_h = B1_h*ninv[0] + B2_h*ninv[1]
            B2n_h = B2n_h / mesure
            hroe  = (w_l["h"]+w_r["h"])/2

            uroe  = un_h
            vroe  = vn_h
            B1roe = B1n_h
            B2roe = B2n_h
            
            w = np.sqrt(pow(B1roe,2) + grav*hroe)
            
            cpsi = max(np.fabs(uroe - w) , np.fabs(uroe + w))

                
            n1     = norm[0]
            n2     = norm[1]
            
            wxl = np.sqrt(pow(B1l,2) + grav*hl)
            wyl = np.sqrt(pow(B2l,2) + grav*hl)
            rx1l = np.fabs(ul - wxl)
            rx2l = np.fabs(ul + wxl)
            ry1l = np.fabs(vl - wyl)
            ry2l = np.fabs(vl + wyl)
            
            wxr = np.sqrt(pow(B1r,2) + grav*hr)
            wyr = np.sqrt(pow(B2r,2) + grav*hr)
            rx1r = np.fabs(ur - wxr)
            rx2r = np.fabs(ur + wxr)
            ry1r = np.fabs(vr - wyr)
            ry2r = np.fabs(vr + wyr)
            
            gammax = max(rx1l, rx2l, rx1r, rx2r,)
            gammay = max(ry1l, ry2l, ry1r, ry2r,)
            

            
            #cpsi   = delta*max(gammax  , gammay )


            
            Unl     = ul*n1  + vl*n2
#            Utl     = vl*n1  - ul*n2    
            Bnl     = B1l*n1 + B2l*n2
#            Btl     = B2l*n1 - B1l*n2
            wl      = np.sqrt(grav * hl + Bnl**2)           

            Unr     = ur*n1  + vr*n2
#            Utr     = vr*n1  - ur*n2    
            Bnr     = B1r*n1 + B2r*n2
 #           Btr     = B2r*n1 - B1r*n2
            wr      = np.sqrt(grav * hr + Bnr**2)      
  
          
#            Unl     = np.sqrt(w_l["hu"]**2 + w_l["hv"]**2)/w_l["h"]
##            Utl     = vl*n1  - ul*n2    
#            Bnl     = np.sqrt(w_l["hB1"]**2 + w_l["hB2"]**2)/w_l["h"]
##            Btl     = B2l*n1 - B1l*n2
#            wl      = np.sqrt(grav * hl + Bnl**2)           
#
#            Unr     = np.sqrt(w_r["hu"]**2 + w_r["hv"]**2)/w_r["h"]
##            Utr     = vr*n1  - ur*n2    
#            Bnr     = np.sqrt(w_r["hB1"]**2 + w_r["hB2"]**2)/w_r["h"]
# #           Btr     = B2r*n1 - B1r*n2
#            wr      = np.sqrt(grav * hr + Bnr**2)                  
    

            lambda1l = np.fabs(Unl - wl)
            lambda2l = np.fabs(Unl - Bnl)
            lambda3l = np.fabs(Unl + Bnl)
            lambda4l = np.fabs(Unl + wl) 
            
            lambda1r = np.fabs(Unr - wr)
            lambda2r = np.fabs(Unr - Bnr)
            lambda3r = np.fabs(Unr + Bnr)
            lambda4r = np.fabs(Unr + wr)
            
            ll = max(lambda1l, lambda2l, lambda3l, lambda4l)
            lr = max(lambda1r, lambda2r, lambda3r, lambda4r)
            
            #cpsi = 1
            lambda_star = max(ll,lr)
            #cpsi = lambda_star
            lambda_star = c_h #max(ll,lr)
            cpsi = c_h
    
            q_l = w_l["hu"] * norm[0] + w_l["hv"] * norm[1]
            m_l = w_l["hB1"] * norm[0] + w_l["hB2"] * norm[1]
            p_l = 0.5 * grav * w_l["h"]* w_l["h"]
            Psl = w_l['hP']
            
            q_r = w_r["hu"] * norm[0] + w_r["hv"] * norm[1]
            m_r = w_r["hB1"] * norm[0] + w_r["hB2"] * norm[1]
            p_r = 0.5 * grav * w_r["h"] * w_r["h"]
            Psr = w_r['hP']

            fleft["h"]  = q_l
            fleft["hu"]  = q_l * w_l["hu"] /w_l["h"]  + p_l*norm[0] - m_l * w_l['hB1'] /w_l["h"]
            fleft["hv"]  = q_l * w_l["hv"] /w_l["h"]  + p_l*norm[1] - m_l * w_l['hB2'] /w_l["h"]
            fleft["hB1"] = (w_l["hv"] * w_l['hB1'] /w_l["h"] - w_l["hu"] * w_l['hB2'] /w_l["h"] )*norm[1] #+ Psl * norm[0]
            fleft["hB2"] = (w_l["hu"] * w_l['hB2'] /w_l["h"] - w_l["hv"] * w_l['hB1'] /w_l["h"])* norm[0] #+ Psl * norm[1]
            fleft["hP"]  = pow(cpsi,2)*m_l
  
            fright["h"] = q_r
            fright["hu"]  = q_r * w_r["hu"] /w_r["h"]  + p_r*norm[0] - m_r * w_r['hB1'] /w_r["h"]
            fright["hv"]  = q_r * w_r["hv"] /w_r["h"]  + p_r*norm[1] - m_r * w_r['hB2'] /w_r["h"]
            fright["hB1"] = (w_r["hv"] * w_r['hB1'] /w_r["h"] - w_r["hu"] * w_r['hB2'] /w_r["h"] )*norm[1] #+ Psr * norm[0]
            fright["hB2"] = (w_r["hu"] * w_r['hB2'] /w_r["h"] - w_r["hv"] * w_r['hB1'] /w_r["h"])* norm[0] #+ Psr * norm[1]
            fright["hP"]  = pow(cpsi,2)*m_r
        
            f_h   = 0.5 * (fleft["h"] + fright["h"])     - 0.5*lambda_star*(w_r["h"]    - w_l["h"])
            f_hu  = 0.5 * (fleft["hu"] + fright["hu"])   - 0.5*lambda_star*(w_r["hu"]   - w_l["hu"])
            f_hv  = 0.5 * (fleft["hv"] + fright["hv"])   - 0.5*lambda_star*(w_r["hv"]   - w_l["hv"])
            f_hB1 = 0.5 * (fleft["hB1"] + fright["hB1"]) - 0.5*lambda_star*(w_r["hB1"]  - w_l["hB1"])
            f_hB2 = 0.5 * (fleft["hB2"] + fright["hB2"]) - 0.5*lambda_star*( w_r["hB2"] - w_l["hB2"])
            f_hP  = 0#0.5 * (fleft["hP"]  + fright["hP"])  - 0.5*cpsi*( w_r["hP"]  - w_l["hP"])

            #######################################################"""
            
            ad=0.0
            uroe  = ad*un_h
            vroe  = ad*vn_h            
            ##############################################################"""
            cpsi = cpsi
            s , pi= ddm.sgn(uroe, 1e-13)
            cc = pow(cpsi,2) - pow(uroe,2)          
            Acpsi = np.zeros((3,3))
            Acpsi[0][0] = 0.0
            Acpsi[0][1] = -vroe*(cpsi - s*uroe)/cc
            Acpsi[0][2] = cpsi
            
            Acpsi[1][0] = 0.0
            Acpsi[1][1] = s 
            Acpsi[1][2] = 0.0  
            
            Acpsi[0][0] = 1/cpsi
            Acpsi[0][1] = vroe*(cpsi*s - uroe)/(cc*cpsi)
            Acpsi[0][2] = 0.0        
   
   
            B1left  = w_l["hB1"]*norm[0] + w_l["hB2"]*norm[1]	
            B2left  = w_l["hB2"]*norm[0] - w_l["hB1"]*norm[1] 
            B1right = w_r["hB1"]*norm[0] + w_r["hB2"]*norm[1]	    
            B2right = w_r["hB2"]*norm[0] - w_r["hB1"]*norm[1] 
               
            w_dif = np.zeros(3)
            w_lrhB1 = (B1left + B1right)/2
            w_lrhB2 = (B2left + B2right)/2
            w_lrhP = (w_l["hP"]  + w_r["hP"])/2
            
            
            w_dif[0] = B1right  - B1left
            w_dif[1] = B2right  - B2left
            w_dif[2] = w_r["hP"] - w_l["hP"]
            
            B1new = 0.
            B2new = 0.
            hPnew = 0.

            for i in range(3):
                B1new += Acpsi[0][i] * w_dif[i]
                B2new += Acpsi[1][i] * w_dif[i]
                hPnew += Acpsi[2][i] * w_dif[i]
            
            
            u_hP  = hPnew/2
            u_hB1 = B1new/2
            u_hB2 = B2new/2
            
            w_lrhP  = w_lrhP  - u_hP
            hw_lrhB1 = w_lrhB1 - u_hB1
            w_lrhB2 = w_lrhB2 - u_hB2
            
            B1new = 0.
            B2new = 0.
            
            B1new = w_lrhB1 * norm[0]  - w_lrhB2 * norm[1]
            B2new = w_lrhB1 * norm[1]  + w_lrhB2 * norm[0] 
            
            #w_lrhB1 = B1new
            #w_lrhB2 = B2new
            p_s = norm[0] * B1new + norm[1] * B2new
            
            w_hP   = w_r["hP"] - w_l["hP"]
            w_hB1 = w_lrhB1/2 - (1/(2*cpsi))*w_dif[2] 
            mhP  = (w_r["hP"] + w_l["hP"])/2 - cpsi*w_dif[0]/2
            
            Flux_B1psi =  mhP * norm[0]
            Flux_B2psi =  mhP * norm[1]
            Flux_hPpsi = cpsi*cpsi*w_hB1
                       
            
            flux["h"]   = f_h * mesure
            flux["hu"]  = f_hu * mesure
            flux["hv"]  = f_hv * mesure
            flux["hP"]  = (f_hP + Flux_hPpsi)  * mesure
            flux["hB1"] = (f_hB1 + Flux_B1psi) * mesure
            flux["hB2"] = (f_hB2 + Flux_B2psi ) * mesure   
            flux["Ch"]  = cpsi         
            flux['Z']   = 0.0


            return flux



#TODO
@njit(fastmath=True)
def explicitscheme_convective(w_c, w_p, w_x, w_y, psi, w_ghost, w_halo, wx_halo, wy_halo, psi_halo, cellidf, faceidc,
                              nodeidc, centerc, cellidc, centerh, mesuref, centerf, normal, halofid,
                              name, ghostcenterf, cellidn, shift, mystruct, order, SCHEME,  grav, c_h, delta, w_cst, w_ghostc ):

    rezidus =np.zeros(len(w_c), dtype=mystruct)
    w_pl  = np.zeros(1, dtype=mystruct)[0]
    w_pr  = np.zeros(1, dtype=mystruct)[0]
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
            w_pl= w_cst[cellidf[i][0]]
            w_pr= w_cst[cellidf[i][1]]
            norm = normal[i]
            mesu = mesuref[i]
            
            if name[i] == 0:
                w_r = w_c[cellidf[i][1]]
                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_l, w_r, norm, mesu, grav,w_pl, w_pr)
                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_l, w_r, norm, mesu, grav,w_pl, w_pr,c_h)
                elif SCHEME=="LF":
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_l, w_r, norm, mesu, grav)
                elif SCHEME=="GLM":
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_l, w_r, norm, mesu, grav, c_h, delta)
                                                      
                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)
                rezidus[cellidf[i][1]] = ddm.add(rezidus[cellidf[i][1]], flx)

            elif name[i] == 5 or name[i] == 6 or name[i] == 7 or name[i] == 8:

                w_r = w_c[cellidf[i][1]]
                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_l, w_r, norm, mesu, grav,w_pl, w_pr)
                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_l, w_r, norm, mesu, grav,w_pl, w_pr,c_h)
                elif SCHEME=="LF":
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_l, w_r, norm, mesu, grav)
                elif SCHEME=="GLM":
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_l, w_r, norm, mesu, grav, c_h, delta)
                                                      
                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)
                
            elif name[i] == 10:
                w_r = w_halo[halofid[i]]
                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_l, w_r, norm, mesu, grav,w_pl, w_pr)
                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_l, w_r, norm, mesu, grav, w_pl, w_pr, c_h)
                elif (SCHEME=="LF"):
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_l, w_r, norm, mesu, grav)                    
                elif (SCHEME=="GLM"):
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_l, w_r, norm, mesu, grav, c_h, delta) 
                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)

            else:
                w_r = w_ghost[i]
                w_pr = w_ghostc[i]
                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_l, w_r, norm, mesu, grav, w_pl, w_pr)                
                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_l, w_r, norm, mesu, grav, w_pl, w_pr, c_h)                 
                elif (SCHEME=="LF"):
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_l, w_r, norm, mesu, grav)
                elif (SCHEME=="GLM"):
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_l, w_r, norm, mesu, grav, c_h, delta)
                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)


    elif order == 2:

        for i in range(nbface):
           
            norm = normal[i]
            mesu = mesuref[i]
            w_pl= w_cst[cellidf[i][0]]
            w_pr= w_cst[cellidf[i][1]]

            if name[i] == 0:
                w_l = w_c[cellidf[i][0]]
                w_r = w_c[cellidf[i][1]]
                w_pl= w_cst[cellidf[i][0]]
                w_pr= w_cst[cellidf[i][1]]
                
                center_left = centerc[cellidf[i][0]]
                center_right = centerc[cellidf[i][1]]

                w_x_left = w_x[cellidf[i][0]]
                w_y_left = w_y[cellidf[i][0]]
                psi_left = psi[cellidf[i][0]]

                w_x_right = w_x[cellidf[i][1]]
                w_y_right = w_y[cellidf[i][1]]
                psi_right = psi[cellidf[i][1]]

                r_l = np.array([centerf[i][0] - center_left[0], centerf[i][1] - center_left[1]])
                w_ln["h"][0]  = w_l["h"]  + psi_left * (w_x_left["h"]  * r_l[0] + w_y_left["h"]  * r_l[1])
                w_ln["hu"][0] = w_l["hu"] + psi_left * (w_x_left["hu"] * r_l[0] + w_y_left["hu"] * r_l[1])
                w_ln["hv"][0] = w_l["hv"] + psi_left * (w_x_left["hv"] * r_l[0] + w_y_left["hv"] * r_l[1])
                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])
                w_ln["hP"][0] = w_l["hP"] + psi_left * (w_x_left["hP"] * r_l[0] + w_y_left["hP"] * r_l[1])
                w_ln["Z"][0]  = w_l["Z"]  + psi_left * (w_x_left["Z"]  * r_l[0] + w_y_left["Z"]  * r_l[1])

                r_r = np.array([centerf[i][0] - center_right[0], centerf[i][1] - center_right[1]])
                w_rn["h"][0]  = w_r["h"]  + psi_right * (w_x_right["h"]  * r_r[0] + w_y_right["h"]  * r_r[1])
                w_rn["hu"][0] = w_r["hu"] + psi_right * (w_x_right["hu"] * r_r[0] + w_y_right["hu"] * r_r[1])
                w_rn["hv"][0] = w_r["hv"] + psi_right * (w_x_right["hv"] * r_r[0] + w_y_right["hv"] * r_r[1])
                w_rn["hB1"][0] = w_r["hB1"] + psi_right * (w_x_right["hB1"] * r_r[0] + w_y_right["hB1"] * r_r[1])
                w_rn["hB2"][0] = w_r["hB2"] + psi_right * (w_x_right["hB2"] * r_r[0] + w_y_right["hB2"] * r_r[1])
                w_rn["hP"][0] = w_r["hP"] + psi_right * (w_x_right["hP"] * r_r[0] + w_y_right["hP"] * r_r[1])
                w_rn["Z"][0]  = w_r["Z"]  + psi_right * (w_x_right["Z"]  * r_r[0] + w_y_right["Z"]  * r_r[1])

                
                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav, w_pl, w_pr)

                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav, w_pl, w_pr, c_h)


                elif (SCHEME=="LF"):
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav)

                elif (SCHEME=="GLM"):
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav, c_h, delta)
                                                        
                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)
                rezidus[cellidf[i][1]] = ddm.add(rezidus[cellidf[i][1]], flx)

    
            elif name[i] == 5 or name[i] == 6 or name[i] == 7 or name[i] == 8:              

                w_l = w_c[cellidf[i][0]]
                w_r = w_c[cellidf[i][1]]
                w_pl= w_cst[cellidf[i][0]]
                w_pr= w_cst[cellidf[i][1]]
                
                center_left = centerc[cellidf[i][0]]
                center_right = centerc[cellidf[i][1]]

                w_x_left = w_x[cellidf[i][0]]
                w_y_left = w_y[cellidf[i][0]]
                psi_left = psi[cellidf[i][0]]

                w_x_right = w_x[cellidf[i][1]]
                w_y_right = w_y[cellidf[i][1]]
                psi_right = psi[cellidf[i][1]]

                r_l = np.array([centerf[i][0] - center_left[0], centerf[i][1] - center_left[1]])
                w_ln["h"][0]  = w_l["h"]  + psi_left * (w_x_left["h"]  * r_l[0] + w_y_left["h"]  * r_l[1])
                w_ln["hu"][0] = w_l["hu"] + psi_left * (w_x_left["hu"] * r_l[0] + w_y_left["hu"] * r_l[1])
                w_ln["hv"][0] = w_l["hv"] + psi_left * (w_x_left["hv"] * r_l[0] + w_y_left["hv"] * r_l[1])
                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])
                w_ln["hP"][0] = w_l["hP"] + psi_left * (w_x_left["hP"] * r_l[0] + w_y_left["hP"] * r_l[1])
                w_ln["Z"][0]  = w_l["Z"]  + psi_left * (w_x_left["Z"]  * r_l[0] + w_y_left["Z"]  * r_l[1])

                r_r = np.array([centerf[i][0] - center_right[0] - shift[cellidf[i][1]][0],
                                centerf[i][1] - center_right[1] - shift[cellidf[i][1]][1]])
                w_rn["h"][0]  = w_r["h"]  + psi_right * (w_x_right["h"]  * r_r[0] + w_y_right["h"]  * r_r[1])
                w_rn["hu"][0] = w_r["hu"] + psi_right * (w_x_right["hu"] * r_r[0] + w_y_right["hu"] * r_r[1])
                w_rn["hv"][0] = w_r["hv"] + psi_right * (w_x_right["hv"] * r_r[0] + w_y_right["hv"] * r_r[1])
                w_rn["hB1"][0] = w_r["hB1"] + psi_right * (w_x_right["hB1"] * r_r[0] + w_y_right["hB1"] * r_r[1])
                w_rn["hB2"][0] = w_r["hB2"] + psi_right * (w_x_right["hB2"] * r_r[0] + w_y_right["hB2"] * r_r[1])
                w_rn["hP"][0] = w_r["hP"] + psi_right * (w_x_right["hP"] * r_r[0] + w_y_right["hP"] * r_r[1])
                w_rn["Z"][0]  = w_r["Z"]  + psi_right * (w_x_right["Z"]  * r_r[0] + w_y_right["Z"]  * r_r[1])

                
                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav, w_pl, w_pr)

                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav, w_pl, w_pr, c_h)

                elif (SCHEME=="LF"):
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav)

                elif (SCHEME=="GLM"):
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav, c_h, delta)
                                        
                
                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)
                #rezidus[cellidf[i][1]] = ddm.add(rezidus[cellidf[i][1]], flx)

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
                w_ln["h"][0]  = w_l["h"]  + psi_left * (w_x_left["h"]  * r_l[0] + w_y_left["h"]  * r_l[1])
                w_ln["hu"][0] = w_l["hu"] + psi_left * (w_x_left["hu"] * r_l[0] + w_y_left["hu"] * r_l[1])
                w_ln["hv"][0] = w_l["hv"] + psi_left * (w_x_left["hv"] * r_l[0] + w_y_left["hv"] * r_l[1])
                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])
                w_ln["hP"][0] = w_l["hP"] + psi_left * (w_x_left["hP"] * r_l[0] + w_y_left["hP"] * r_l[1])
                w_ln["Z"][0]  = w_l["Z"]  + psi_left * (w_x_left["Z"]  * r_l[0] + w_y_left["Z"]  * r_l[1])

                r_r = np.array([centerf[i][0] - centerh[halofid[i]][0],
                                centerf[i][1] - centerh[halofid[i]][1]])
                w_rn["h"][0]  = w_r["h"]  + psi_right * (w_x_right["h"]  * r_r[0] + w_y_right["h"]  * r_r[1])
                w_rn["hu"][0] = w_r["hu"] + psi_right * (w_x_right["hu"] * r_r[0] + w_y_right["hu"] * r_r[1])
                w_rn["hv"][0] = w_r["hv"] + psi_right * (w_x_right["hv"] * r_r[0] + w_y_right["hv"] * r_r[1])
                w_rn["hB1"][0] = w_r["hB1"] + psi_right * (w_x_right["hB1"] * r_r[0] + w_y_right["hB1"] * r_r[1])
                w_rn["hB2"][0] = w_r["hB2"] + psi_right * (w_x_right["hB2"] * r_r[0] + w_y_right["hB2"] * r_r[1])
                w_rn["hP"][0] = w_r["hP"] + psi_right * (w_x_right["hP"] * r_r[0] + w_y_right["hP"] * r_r[1])
                w_rn["Z"][0]  = w_r["Z"]  + psi_right * (w_x_right["Z"]  * r_r[0] + w_y_right["Z"]  * r_r[1])

                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav,w_pl, w_pr)

                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav,w_pl, w_pr, c_h)
                
                elif (SCHEME=="LF"):
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav)
                elif (SCHEME=="GLM"):
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_ln[0], w_rn[0], norm, mesu, grav, c_h, delta)

                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)

            else:
                w_l = w_c[cellidf[i][0]]
                w_r = w_ghost[i]
                w_pr = w_ghostc[i]

                center_left = centerc[cellidf[i][0]]
                center_right = ghostcenterf[i]

                w_x_left = w_x[cellidf[i][0]]
                w_y_left = w_y[cellidf[i][0]]
                psi_left = psi[cellidf[i][0]]

                r_l = np.array([centerf[i][0] - center_left[0], centerf[i][1] - center_left[1]])
                w_ln["h"][0]  = w_l["h"]  + psi_left * (w_x_left["h"]  * r_l[0] + w_y_left["h"]  * r_l[1])
                w_ln["hu"][0] = w_l["hu"] + psi_left * (w_x_left["hu"] * r_l[0] + w_y_left["hu"] * r_l[1])
                w_ln["hv"][0] = w_l["hv"] + psi_left * (w_x_left["hv"] * r_l[0] + w_y_left["hv"] * r_l[1])
                w_ln["hB1"][0] = w_l["hB1"] + psi_left * (w_x_left["hB1"] * r_l[0] + w_y_left["hB1"] * r_l[1])
                w_ln["hB2"][0] = w_l["hB2"] + psi_left * (w_x_left["hB2"] * r_l[0] + w_y_left["hB2"] * r_l[1])
                w_ln["hP"][0] = w_l["hP"] + psi_left * (w_x_left["hP"] * r_l[0] + w_y_left["hP"] * r_l[1])
                w_ln["Z"][0]  = w_l["Z"]  + psi_left * (w_x_left["Z"]  * r_l[0] + w_y_left["Z"]  * r_l[1])

                if SCHEME=="SRNH":
                    flx = compute_flux_shallow_srnh_mhd(flx, fleft, fright, w_ln[0], w_r, norm, mesu, grav, w_pl, w_pr)                

                elif SCHEME=="SRNH-GLM":
                    flx = compute_flux_shallow_srnh_mhd_GLM(flx, fleft, fright, w_ln[0], w_r, norm, mesu, grav, w_pl, w_pr,c_h)                


                elif (SCHEME=="LF"):
                    flx = compute_flux_swmhd_LF(flx, fleft, fright, w_ln[0], w_r, norm, mesu, grav)

                elif (SCHEME=="GLM"):
                    flx = compute_flux_swmhd_LF_GLM_correction(flx, fleft, fright, w_ln[0], w_r, norm, mesu, grav, c_h,delta)
            
                                
                rezidus[cellidf[i][0]] = ddm.minus(rezidus[cellidf[i][0]], flx)


                
    return rezidus


@njit(fastmath=True)
def term_coriolis(w_c, f_c, mystruct):
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
    
    
    
    

 #TODO source
@njit(fastmath=True)
def term_source_srnh_SW(src, w_c, w_ghost, w_halo, w_x, w_y, psi, wx_halo, wy_halo, psi_halo, 
                        nodeidc, faceidc, cellidc,  cellidf, centerc, normalc, namef, centerf, centerh,
                        vertexn, halofid, grav, order):

 
    nbelement = len(w_c)
    hi_p  =  np.zeros(3)
    zi_p  =  np.zeros(3)
    B1i_p =  np.zeros(3)
    B2i_p =  np.zeros(3)
#    r_l = zeros(2)
#    r_r = zeros(2)

    zv =  np.zeros(3)
    
    mata =  np.zeros(3)
    matb =  np.zeros(3)
    
    ns =  np.zeros((3, 3))
    ss =  np.zeros((3, 3))
    s_1 = np.zeros(3)
    s_2 = np.zeros(3)
    s_3 = np.zeros(3)
    b = np.zeros(3)


    for i in range(nbelement):
        
        G = centerc[i]
        c_1 = 0.
        c_2 = 0.

        for j in range(3):
            f = faceidc[i][j]
            ss[j] = normalc[i][j]
            
            if namef[f] == 10 :
                
                h_1p  = w_c.h[i]
                z_1p  = w_c.Z[i]
                B1_1p = w_c.hB1[i]
                B2_1p = w_c.hB2[i]
                
                h_p1 = w_halo.h[halofid[f]]
                z_p1 = w_halo.Z[halofid[f]]
                B1_p1 = w_halo.hB1[halofid[f]]
                B2_p1 = w_halo.hB2[halofid[f]]
                
            elif namef[f] == 0:
            
                h_1p = w_c.h[i]
                z_1p = w_c.Z[i]
                B1_1p = w_c.hB1[i]
                B2_1p = w_c.hB2[i]
                                
                h_p1  = w_c.h[cellidc[i][j]]
                z_p1  = w_c.Z[cellidc[i][j]]
                B1_p1 = w_c.hB1[cellidc[i][j]]
                B2_p1 = w_c.hB2[cellidc[i][j]]

            elif namef[i] == 5 or namef[i] == 6 or namef[i] == 7 or namef[i] == 8 :
                h_1p = w_c.h[i]
                z_1p = w_c.Z[i]
                B1_1p = w_c.hB1[i]
                B2_1p = w_c.hB2[i]
                                
                h_p1  = w_c.h[cellidc[i][j]]
                z_p1  = w_c.Z[cellidc[i][j]]
                B1_p1 = w_c.hB1[cellidc[i][j]]
                B2_p1 = w_c.hB2[cellidc[i][j]]             
            else:
                h_1p = w_c.h[i]
                z_1p = w_c.Z[i]
                B1_1p = w_c.hB1[i]
                B2_1p = w_c.hB2[i]
                
                h_p1 = w_ghost.h[f]
                z_p1 = w_ghost.Z[f]
                B1_p1 = w_ghost.hB1[f]
                B2_p1 = w_ghost.hB2[f]
                
            zv[j] = z_p1
            mata[j] = h_p1*ss[j][0]
            matb[j] = h_p1*ss[j][1]
            
            hij  = 0.5*(h_1p + h_p1)
#            B1ij = 0.5*(B1_1p/h_1p + B1_p1/h_p1)
#            B2ij = 0.5*(B2_1p/h_1p + B2_p1/h_p1)
            B1ij = 0.5*(B1_1p + B1_p1)/hij
            B2ij = 0.5*(B2_1p + B2_p1)/hij
            
            c_1 = c_1 + (0.5*(h_1p + h_p1)*0.5*(h_1p + h_p1))  *ss[j][0]
            c_2 = c_2 + (0.5*(h_1p + h_p1)*0.5*(h_1p + h_p1))  *ss[j][1]
            
            #c_1 = c_1 + (pow(hij,2) - 2*hij*pow(B1ij,2))*ss[j][0] - 2*hij*B1ij*B2ij*ss[j][1] 
            #c_2 = c_2 + (pow(hij,2) - 2*hij*pow(B2ij,2))*ss[j][1] - 2*hij*B1ij*B2ij*ss[j][0]            
            
            hi_p[j] = h_1p
            zi_p[j] = z_1p
            

        c_3 = 3.0 * h_1p
            
        delta = (mata[1]*matb[2]-mata[2]*matb[1]) - (mata[0]*matb[2]-matb[0]*mata[2]) + (mata[0]*matb[1]-matb[0]*mata[1])

        deltax = c_3*(mata[1]*matb[2]-mata[2]*matb[1]) - (c_1*matb[2]-c_2*mata[2]) + (c_1*matb[1]-c_2*mata[1])

        deltay = (c_1*matb[2]-c_2*mata[2]) - c_3*(mata[0]*matb[2]-matb[0]*mata[2]) + (mata[0]*c_2-matb[0]*c_1)

        deltaz = (mata[1]*c_2-matb[1]*c_1) - (mata[0]*c_2-matb[0]*c_1) + c_3*(mata[0]*matb[1]-matb[0]*mata[1])
        
        h_1 = deltax/delta
        h_2 = deltay/delta
        h_3 = deltaz/delta
            
        z_1 = zi_p[0] + hi_p[0] - h_1
        z_2 = zi_p[1] + hi_p[1] - h_2
        z_3 = zi_p[2] + hi_p[2] - h_3

        b =  np.array([vertexn[nodeidc[i][1]][0], vertexn[nodeidc[i][1]][1], 0.])

        ns[0] = np.array([(G[1]-b[1]), -(G[0]-b[0]), 0.])
        ns[1] = ns[0] - ss[1]  # N23
                                                                                                                                                                      
        ns[2] = ns[0] + ss[0]  #  N31    
        
        s_1 = 0.5*h_1*(zv[0]*ss[0] + z_2*ns[0] + z_3*(-1)*ns[2])
        s_2 = 0.5*h_2*(zv[1]*ss[1] + z_1*(-1)*ns[0] + z_3*ns[1])
        s_3 = 0.5*h_3*(zv[2]*ss[2] + z_1*ns[2] + z_2*(-1)*ns[1])
        
        #TODO
        src[i]["h"] = 0
        src[i]["hu"] = -grav*( (s_1[0] + s_2[0] + s_3[0]) )
        src[i]["hv"] = -grav*( (s_1[1] + s_2[1] + s_3[1]) )
        src[i]["hB1"] = 0.
        src[i]["hB2"] = 0.
        src[i]["hP"] = 0.
        src[i]["Ch"] = 0.
        src[i]["Z"] = 0.



