from pyccel.decorators import stack_array, inline


#@stack_array('norm') 
#def explicitscheme_stabilizater(ux_face:'float[:]',  uy_face:'float[:]', vx_face:'float[:]',  vy_face:'float[:]',
#                               cellidf:'int[:,:]', normalf:'float[:,:]', namef:'int[:]', 
#                               dissip_u:'float[:]', dissip_v:'float[:]', h_c:'float[:]', epsilon:'float'):
#    
#    from numpy import zeros
#    
#    nbface = len(cellidf)
#    norm = zeros(3)
#    dissip_u[:] = 0.
#    dissip_v[:] = 0.
#
#    for i in range(nbface):
#        
#        norm[:] = normalf[i][:]
#        q = epsilon*h_c[i] * (ux_face[i] * norm[0] + uy_face[i] * norm[1])
#        p = epsilon*h_c[i] * (vx_face[i] * norm[0] + vy_face[i] * norm[1])
#
#        flux_u = q
#        flux_v = p
#
#        if namef[i] == 0:
#
#            dissip_u[cellidf[i][0]] += flux_u
#            dissip_u[cellidf[i][1]] -= flux_u
#            dissip_v[cellidf[i][0]] += flux_v
#            dissip_v[cellidf[i][1]] -= flux_v
#        else:
#            dissip_u[cellidf[i][0]] += flux_u
#            dissip_v[cellidf[i][0]] += flux_v


def Total_Energy(h_c:'float[:]', hu_c:'float[:]', hv_c:'float[:]', hB1_c:'float[:]', hB2_c:'float[:]', Z_c:'float[:]', grav:'float', volumec: 'float[:]'):

    from math import pow
    from numpy import zeros
    nbelement = len(h_c)
    num_t = 0. ; num_c = 0. ; num_m = 0. ; num_p = 0.
    #from numpy import zeros
    Et = zeros(nbelement) # Total_Energy
    Ec = zeros(nbelement) # Kinetic Energy
    Ep = zeros(nbelement) # Potential Energy
    Em = zeros(nbelement) # Magnetic energy
    
    for i in range(nbelement):
        hc  = h_c[i]
        uc  = hu_c[i]/h_c[i]
        vc  = hv_c[i]/h_c[i]
        B1c = hB1_c[i]/h_c[i]
        B2c = hB2_c[i]/h_c[i]
        bc  = Z_c[i]
        Ec[i] = 0.5*hc*(pow(uc , 2)  + pow(vc , 2))
        Em[i] = 0.5*hc*(pow(B1c , 2) + pow(B2c , 2))
        Ep[i] = 0.5*grav*pow(hc , 2) + grav*bc*hc
        Et[i] = Ec[i] + Em[i] + Ep[i] 
        num_t += volumec[i]*Et[i] 
        num_c += volumec[i]*Ec[i] 
        num_m += volumec[i]*Em[i] 
        num_p += volumec[i]*Ep[i] 

    numt = num_t ; numc = num_c ; numm = num_m ; nump = num_p 
    return numt , numc , numm , nump   

@stack_array('norm') 
def explicitscheme_stabilizater(wx_face:'float[:]',  wy_face:'float[:]', cellidf:'int[:,:]', normalf:'float[:,:]', 
                                namef:'int[:]', vepsilon:'float', dissip_w:'float[:]', h_c:'float[:]'):
    
    from numpy import zeros
    
    nbface = len(cellidf)
    norm = zeros(3)
    dissip_w[:] = 0.

    for i in range(nbface):
        
        norm[:] = normalf[i][:]
        q = vepsilon*h_c[i] * (wx_face[i] * norm[0] + wy_face[i] * norm[1])

        flux_w = q

        if namef[i] == 0:

            dissip_w[cellidf[i][0]] += flux_w
            dissip_w[cellidf[i][1]] -= flux_w

        else:
            dissip_w[cellidf[i][0]] += flux_w

#renvoie min
@stack_array('norm')          
def time_step_SWMHD(h_c:'float[:]', hu_c:'float[:]', hv_c:'float[:]', hB1_c:'float[:]', hB2_c:'float[:]',
                  cfl:'float', normal:'float[:,:]', mesure:'float[:]', volume:'float[:]', 
                  faceid:'int[:,:]'):
   
    from numpy import sqrt, fabs
    grav = 1.0
    nbelement =  len(faceid)
    u_n = 0.
    B_n = 0.
    dt = 1e6
 
    for i in range(nbelement):
        lam = 0.
        for j in range(3):
            u_n = (fabs(hu_c[i]/h_c[i]*normal[faceid[i][j]][0] + hv_c[i]/h_c[i]*normal[faceid[i][j]][1]))/mesure[faceid[i][j]] 
            B_n = (fabs(hB1_c[i]/h_c[i]*normal[faceid[i][j]][0] + hB2_c[i]/h_c[i]*normal[faceid[i][j]][1]))/mesure[faceid[i][j]] 
            wb  = sqrt(B_n**2 + grav*h_c[i])
            lam1 = fabs(u_n - wb)
            lam2 = fabs(u_n - B_n)
            lam3 = fabs(u_n + B_n)
            lam4 = fabs(u_n + wb)
            #lam_convect = u_n + velson #max(lam1, lam2, lam3, lam4)
            lam_convect = max(lam1, lam2, lam3, lam4) 
            lam += lam_convect * mesure[faceid[i][j]]
            #lam += lam_convect * mesure[faceid[i][j]]
                 
       # dt_c[i]  = cfl * volume[i]/lam

        dt  = min(dt, cfl * volume[i]/lam)
        
    return dt
def update_SWMHD(h_c:'float[:]', hu_c:'float[:]', hv_c:'float[:]', hB1_c:'float[:]', hB2_c:'float[:]', PSI_c:'float[:]', Z_c:'float[:]', 
              rez_h:'float[:]', rez_hu:'float[:]', rez_hv:'float[:]', rez_hB1:'float[:]', rez_hB2:'float[:]', rez_PSI:'float[:]', rez_Z:'float[:]',
              src_h:'float[:]', src_hu:'float[:]', src_hv:'float[:]', src_hB1:'float[:]', src_hB2:'float[:]', src_PSI:'float[:]', src_Z:'float[:]',
              dtime:'float', vol:'float[:]', GLM:'int', cpsi: 'float'):

    
    from numpy import exp
    h_c[:]   += dtime  * (rez_h[:]    + src_h[:] )/vol[:]
    hu_c[:]  += dtime  * (rez_hu[:]  + src_hu[:] )/vol[:] 
    hv_c[:]  += dtime  * (rez_hv[:]  + src_hv[:] )/vol[:] 
    hB1_c[:] += dtime  * (rez_hB1[:]  + src_hB1[:] )/vol[:]
    hB2_c[:] += dtime  * (rez_hB2[:]  + src_hB2[:] )/vol[:]
    PSI_c[:] += dtime  * (rez_PSI[:]  + src_PSI[:] )/vol[:]
    Z_c[:]   += dtime  * (rez_Z[:]    + src_Z[:] )/vol[:]

    if GLM==10:
         cr = 0.01
         PSI_c[:] = exp(-dtime*(cpsi/cr))*PSI_c[:]        
  
  #renvoi max  
@stack_array('norm')          
def cpsi_global(h_c:'float[:]', hu_c:'float[:]', hv_c:'float[:]', hB1_c:'float[:]', hB2_c:'float[:]',
                  cfl:'float', normal:'float[:,:]', mesure:'float[:]', volume:'float[:]', 
                  faceid:'int[:,:]'):
   
    from numpy import sqrt, fabs
    grav = 1
    nbelement =  len(faceid)
    u_n = 0.
    B_n = 0. 
    cpsiglobal = 0.
    
    
    for i in range(nbelement):
        lam = 0.
        for j in range(3):
            u_n = (fabs(hu_c[i]/h_c[i]*normal[faceid[i][j]][0] + hv_c[i]/h_c[i]*normal[faceid[i][j]][1]))/mesure[faceid[i][j]] 
            B_n = (fabs(hB1_c[i]/h_c[i]*normal[faceid[i][j]][0] + hB2_c[i]/h_c[i]*normal[faceid[i][j]][1]))/mesure[faceid[i][j]] 
            wb  = sqrt(B_n**2 + grav*h_c[i])
            lam1 = fabs(u_n - wb)
            lam2 = fabs(u_n - B_n)
            lam3 = fabs(u_n + B_n)
            lam4 = fabs(u_n + wb)
            lam_convect = max(lam1, lam2, lam3, lam4)
            lam += lam_convect * mesure[faceid[i][j]]
        cpsiglobal  = max(cpsiglobal, lam_convect)
        
    return cpsiglobal


def initialisation_SWMHD(h:'float[:]', hu:'float[:]', hv:'float[:]', hB1:'float[:]', hB2:'float[:]', PSI:'float[:]', Z:'float[:]',
                      center:'float[:,:]', choix:'int', k1:'float', k2:'float', eps:'float', tol:'float'):
   
    from numpy import  sqrt, pi, sin, exp, cos, arccos, fabs
    from math import pow
     #fabs, pi, cos, sin, exp
    
    nbelements = len(center)
    
    

    if choix == 1:  
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            hint = 1./4
            uint = 1 + 0.5*sin(pi*ycent) + hint*cos(pi*xcent)
            vint = 1 + hint*sin(pi*xcent) + 0.5*cos(pi*ycent)


            h[i]  = hint
            hu[i] =  hint*uint 
            hv[i] =  hint*vint 
            hB1[i] = hint*(0.5)
            hB2[i] = hint*uint 





    if choix == 7: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            AA = sin(2*pi*xcent-2*pi*ycent)

            uu  = 2.
            bb  = 5.
            h0  = 100.
            u0  = (uu + AA)/h0
            v0  = (uu + AA)/h0
            B01 = (bb + 2*AA)/h0
            B02 = (bb + 2*AA)/h0



            PSI[i] = 0. 
            h[i]  = 100
            hu[i] = uu + AA
            hv[i] = uu + AA
            hB1[i] = bb + 2*AA
            hB2[i] = bb + 2*AA  
            PSI[i] = 0.
            
            grav = 1.0            
            k  = sqrt(k1**2 + k2**2)
            #B0 = sqrt(B01**2 + B02**2)  
            c0 = sqrt(h0*grav)
            #Produit scalaire de (k1,k2) et (B01, B02)
            ss = k1*B01 + k2*B02
            a1 = sqrt(ss**2 + (k*c0)**2)
            b1 = - 0.5*eps*k*k


            h_hot  =  h0*k*k*cos(k1*xcent + k2*ycent)
            u_hot  =  k1*(a1*cos(k1*xcent + k2*ycent) - b1*sin(k1*xcent + k2*ycent))
            v_hot  =  k2*(a1*cos(k1*xcent + k2*ycent) - b1*sin(k1*xcent + k2*ycent))
            B1_hot = -k1*ss*cos(k1*xcent + k2*ycent)
            B2_hot = -k2*ss*cos(k1*xcent + k2*ycent)

            
            hint  = h0  + tol*h_hot #*cos(k1*xcent + k2*ycent)
            uint  = u0  + tol*u_hot #* cos(k1*xcent + k2*ycent)
            vint  = v0  + tol*v_hot #* cos(k1*xcent + k2*ycent)
            B1int = B01 + tol*B1_hot#* cos(k1*xcent + k2*ycent)
            B2int = B02 + tol*B2_hot#* cos(k1*xcent + k2*ycent)
      
            h[i]  = hint
            hu[i] =  hint*uint 
            hv[i] =  hint*vint 
            hB1[i] = hint*B1int 
            hB2[i] = hint*B2int 


    if choix == -3: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            AA = 0.0*sin(2*pi*xcent-2*pi*ycent)

            uu  = 0.
            bb  = 3.
            h0  = 2.
            u0  = (uu + AA)/h0
            v0  = (uu + AA)/h0
            B01 = (bb + 2*AA)/h0
            B02 = (bb + 2*AA)/h0



            PSI[i] = 0. 
            h[i]  = 100
            hu[i] = uu + AA
            hv[i] = uu + AA
            hB1[i] = bb + 2*AA
            hB2[i] = bb + 2*AA  
            PSI[i] = 0.
            
            grav = 1.0            
            k  = sqrt(k1**2 + k2**2)
            #B0 = sqrt(B01**2 + B02**2)  
            c0 = sqrt(h0*grav)
            #Produit scalaire de (k1,k2) et (B01, B02)
            ss = k1*B01 + k2*B02
            a1 = sqrt(ss**2 + (k*c0)**2)
            b1 = - 0.5*eps*k*k


            h_hot  =  h0*k*k*cos(k1*xcent + k2*ycent)
            u_hot  =  k1*(a1*cos(k1*xcent + k2*ycent) - b1*sin(k1*xcent + k2*ycent))
            v_hot  =  k2*(a1*cos(k1*xcent + k2*ycent) - b1*sin(k1*xcent + k2*ycent))
            B1_hot = -k1*ss*cos(k1*xcent + k2*ycent)
            B2_hot = -k2*ss*cos(k1*xcent + k2*ycent)

            
            hint  = h0  + tol*h_hot #*cos(k1*xcent + k2*ycent)
            uint  = u0  + tol*u_hot #* cos(k1*xcent + k2*ycent)
            vint  = v0  + tol*v_hot #* cos(k1*xcent + k2*ycent)
            B1int = B01 + tol*B1_hot#* cos(k1*xcent + k2*ycent)
            B2int = B02 + tol*B2_hot#* cos(k1*xcent + k2*ycent)
      
            h[i]  = hint
            hu[i] =  hint*uint 
            hv[i] =  hint*vint 
            hB1[i] = hint*B1int 
            hB2[i] = hint*B2int 






    if choix == -2: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            AA = sin(2*pi*xcent-2*pi*ycent)

            uu  = 2.
            bb  = 5.
            h0  = 100.
            u0  = (uu + AA)/h0
            v0  = (uu + AA)/h0
            B01 = (bb + 2*AA)/h0
            B02 = (bb + 2*AA)/h0



            PSI[i] = 0. 
            h[i]  = 100
            hu[i] = uu + AA
            hv[i] = uu + AA
            hB1[i] = bb + 2*AA
            hB2[i] = bb + 2*AA  
            PSI[i] = 0.
            
            grav = 1.0            
            k  = sqrt(k1**2 + k2**2)
            #B0 = sqrt(B01**2 + B02**2)  
            c0 = sqrt(h0*grav)
            #Produit scalaire de (k1,k2) et (B01, B02)
            ss = k1*B01 + k2*B02
            b1 = - 0.5*eps*k*k


            h_hot  =  0.
            u_hot  =  k2*ss*cos(k1*xcent + k2*ycent) 
            v_hot  = -k1*ss*cos(k1*xcent + k2*ycent)
            B1_hot = -k2*(ss*cos(k1*xcent + k2*ycent) + b1*sin(k1*xcent + k2*ycent))
            B2_hot = -k1*(ss*cos(k1*xcent + k2*ycent) + b1*sin(k1*xcent + k2*ycent))

            
            hint  = h0  + tol*h_hot #*cos(k1*xcent + k2*ycent)
            uint  = u0  + tol*u_hot #* cos(k1*xcent + k2*ycent)
            vint  = v0  + tol*v_hot #* cos(k1*xcent + k2*ycent)
            B1int = B01 + tol*B1_hot#* cos(k1*xcent + k2*ycent)
            B2int = B02 + tol*B2_hot#* cos(k1*xcent + k2*ycent)
      
            h[i]  = hint
            hu[i] =  hint*uint 
            hv[i] =  hint*vint 
            hB1[i] = hint*B1int 
            hB2[i] = hint*B2int 

           









    elif choix == 40: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            g = 1.0  ;   umax = 0.2  ;   Bmax = 0.1  ;   hmax = 1.0
            
            rcent = sqrt(xcent**2 + ycent**2)            
            ee    = exp(1-rcent**2)
            e1    = exp(0.5*(1-rcent**2))
            hin   = hmax - (1./(2.0*g))*(umax**2 - Bmax**2)*ee
            uin   = 1.0 - umax*e1*ycent
            vin   = 1.0 + umax*e1*xcent
            B1in  = -Bmax*e1*ycent
            B2in  = Bmax*e1*xcent
            
            h[i]  = hin
            hu[i] = hin*uin
            hv[i] = hin*vin
            hB1[i] = hin*B1in
            hB2[i] = hin*B2in  
            PSI[i] = 0.   




    elif choix == 60: 
    
        for i in range(nbelements):
            xcent = center[i][0]
            ycent = center[i][1]
            
            AA = sin(2*pi*xcent-2*pi*ycent)

            uu = 2.
            bb = 5.
            h[i]  = 100
            hu[i] = uu + AA
            hv[i] = uu + AA
            hB1[i] = bb + 2*AA
            hB2[i] = bb + 2*AA  
            PSI[i] = 0.   




 


        
def term_source_srnh_SWMHD(src_h:'float[:]', src_hu:'float[:]', src_hv:'float[:]', src_hB1:'float[:]', src_hB2:'float[:]', src_PSI:'float[:]', src_Z:'float[:]', 
                        h_c:'float[:]', hu_c:'float[:]', hv_c:'float[:]', Z_c:'float[:]', 
                        h_ghost:'float[:]', hu_ghost:'float[:]', hv_ghost:'float[:]', Z_ghost:'float[:]',
                        h_halo:'float[:]', hu_halo:'float[:]', hv_halo:'float[:]', Z_halo:'float[:]',
                        h_x:'float[:]', h_y:'float[:]', psi:'float[:]',
                        hx_halo:'float[:]', hy_halo:'float[:]', psi_halo:'float[:]', 
                        nodeidc:'int[:,:]', faceidc:'int[:,:]', cellidc:'int[:,:]',  cellidf:'int[:,:]',
                        centerc:'float[:,:]', normalc:'float[:,:,:]', 
                        namef:'int[:]', centerf:'float[:,:]', centerh:'float[:,:]',
                        vertexn:'float[:,:]', halofid:'int[:]', order:'int'):
    
    from numpy import zeros, array

    grav = 9.81
    nbelement = len(h_c)
    hi_p =  zeros(3)
    zi_p =  zeros(3)

    zv =  zeros(3)
    
    mata =  zeros(3)
    matb =  zeros(3)
    
    ns  = zeros((3, 3))
    ss  = zeros((3, 3))
    s_1 = zeros(3)
    s_2 = zeros(3)
    s_3 = zeros(3)
    b   = zeros(3)
    G   = zeros(3)

    for i in range(nbelement):
        
        G[:] = centerc[i]
        c_1 = 0.
        c_2 = 0.

        for j in range(3):
            f = faceidc[i][j]
            ss[j] = normalc[i][j]
            
            if namef[f] == 10 :
                
                h_1p = h_c[i]
                z_1p = Z_c[i]
                
                h_p1 = h_halo[halofid[f]]
                z_p1 = Z_halo[halofid[f]]
                
            elif namef[f] == 0:
                
                h_1p = h_c[i]
                z_1p = Z_c[i]
                
                h_p1 = h_c[cellidc[i][j]]
                z_p1 = Z_c[cellidc[i][j]]
                
            else:
                h_1p = h_c[i]
                z_1p = Z_c[i]
                
                h_p1 = h_ghost[f]
                z_p1 = Z_ghost[f]
                
            zv[j] = z_p1
            mata[j] = h_p1*ss[j][0]
            matb[j] = h_p1*ss[j][1]

            
            c_1 = c_1 + (0.5*(h_1p + h_p1)*0.5*(h_1p + h_p1))  *ss[j][0]
            c_2 = c_2 + (0.5*(h_1p + h_p1)*0.5*(h_1p + h_p1))  *ss[j][1]
            
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

        b[:] =  vertexn[nodeidc[i][1]][0:3]

        ns[0] = array([(G[1]-b[1]), -(G[0]-b[0]), 0.])
        ns[1] = ns[0] - ss[1]  #  N23                                                                                                                                                                       
        ns[2] = ns[0] + ss[0]  #  N31    
        
        s_1 = 0.5*h_1*(zv[0]*ss[0] + z_2*ns[0] + z_3*(-1)*ns[2])
        s_2 = 0.5*h_2*(zv[1]*ss[1] + z_1*(-1)*ns[0] + z_3*ns[1])
        s_3 = 0.5*h_3*(zv[2]*ss[2] + z_1*ns[2] + z_2*(-1)*ns[1])
        
        #TODO
        src_h[i]   = 0.
        src_hu[i]  = -grav * (s_1[0] + s_2[0] + s_3[0])
        src_hv[i]  = -grav * (s_1[1] + s_2[1] + s_3[1]) 
        src_hB1[i]  = 0.
        src_hB2[i]  = 0.
        src_PSI[i]  = 0.
        src_Z[i]   = 0.


@inline
def LF_scheme_MHD(hu_l:'float', hu_r:'float', hv_l:'float', hv_r:'float', h_l:'float', h_r:'float', hB1_l:'float', hB1_r:'float', hB2_l:'float',hB2_r:'float',
                  Z_l:'float', Z_r:'float', normal:'float[:]', mesure:'float', grav:'float', flux:'float[:]'):
    
    
    
    from numpy import zeros, sqrt, fabs
#@inline
#def compute_flux_swmhd_LF(flux, fleft, fright, w_l, w_r, normal, mesure, grav):
        #print("Im Roe scheme")
    norm = normal/mesure
    ninv =zeros(2)
    ninv[0] = -1*norm[1]
    ninv[1] = norm[0]

        
    hl  = h_l   
    ul  = hu_l/ h_l
    vl  = hv_l/ h_l
    B1l = hB1_l/ h_l
    B2l = hB2_l/ h_l
    hr  =  h_r   
    ur  = hu_r/ h_r
    vr  = hv_r/ h_r
    B1r = hB1_r/ h_r
    B2r = hB2_r/ h_r
    
    u_h = (hu_l/ h_l* sqrt(h_l)
	          + hu_r/ h_r* sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))
    v_h = (hv_l/ h_l* sqrt(h_l)
	    	   + hv_r/ h_r* sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))
    B1_h = (hB1_l/ h_l* sqrt(h_l)
	       + hB1_r/ h_r* sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))
	    
    B2_h = (hB2_l/ h_l* sqrt(h_l)
		   + hB2_r/ h_r* sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))  
    un_h = u_h*normal[0] + v_h*normal[1] 
    un_h = un_h / mesure
    vn_h = u_h*ninv[0] + v_h*ninv[1]
    vn_h = vn_h / mesure
    B1n_h = B1_h*normal[0] + B2_h*normal[1]
    B1n_h = B1n_h / mesure
    B2n_h = B1_h*ninv[0] + B2_h*ninv[1]
    B2n_h = B2n_h / mesure
    hroe  = (h_l+h_r)/2
    uroe  = un_h
    vroe  = vn_h
    B1roe = B1n_h
    B2roe = B2n_h
    n1     = norm[0]
    n2     = norm[1]
    
    Unl     = ul*n1  + vl*n2   
    Bnl     = B1l*n1 + B2l*n2
    wl      = sqrt(grav * hl + Bnl**2)
    Unr     = ur*n1  + vr*n2
    Bnr     = B1r*n1 + B2r*n2
    wr      = sqrt(grav * hr + Bnr**2)
        
#    Utr     = vr*n1  - ur*n2 
#    Btr     = B2r*n1 - B1r*n2   
#    Utl     = vl*n1  - ul*n2 
#    Btl     = B2l*n1 - B1l*n2
             
    

    lambda1l = fabs(Unl - wl)
    lambda2l = fabs(Unl - Bnl)
    lambda3l = fabs(Unl + Bnl)
    lambda4l = fabs(Unl + wl) 
    
    lambda1r = fabs(Unr - wr)
    lambda2r = fabs(Unr - Bnr)
    lambda3r = fabs(Unr + Bnr)
    lambda4r = fabs(Unr + wr)
    
    ll = max(lambda1l, lambda2l, lambda3l, lambda4l)
    lr = max(lambda1r, lambda2r, lambda3r, lambda4r)
    
    lambda_star = max(ll,lr)
    
    q_l = hu_l* norm[0] + hv_l* norm[1]
    m_l = hB1_l* norm[0] + hB2_l* norm[1]
    p_l = 0.5 * grav * h_l* h_l
    
    q_r = hu_r* norm[0] + hv_r* norm[1]
    m_r = hB1_r* norm[0] + hB2_r* norm[1]
    p_r = 0.5 * grav * h_r* h_r

    fleft_h  = q_l
    fleft_hu  = q_l * hu_l/h_l + p_l*norm[0] - m_l * hB1_l /h_l
    fleft_hv  = q_l * hv_l/h_l + p_l*norm[1] - m_l * hB2_l /h_l
    fleft_hB1 = (hv_l* hB1_l /h_l- hu_l* hB2_l /h_l)*norm[1]
    fleft_hB2 = (hu_l* hB2_l /h_l- hv_l* hB1_l /h_l)* norm[0]
  
    fright_h = q_r
    fright_hu  = q_r * hu_r/h_r + p_r*norm[0] - m_r * hB1_r /h_r
    fright_hv  = q_r * hv_r/h_r + p_r*norm[1] - m_r * hB2_r /h_r
    fright_hB1 = (hv_r* hB1_r /h_r- hu_r* hB2_r /h_r)*norm[1]
    fright_hB2 = (hu_r* hB2_r /h_r- hv_r* hB1_r /h_r)* norm[0]
        
    f_h = 0.5 * (fleft_h + fright_h)       - 0.5*lambda_star*(h_r  - h_l)
    f_hu = 0.5 * (fleft_hu + fright_hu)    - 0.5*lambda_star*(hu_r - hu_l)
    f_hv = 0.5 * (fleft_hv + fright_hv)    - 0.5*lambda_star*(hv_r - hv_l)
    f_hB1 = 0.5 * (fleft_hB1 + fright_hB1) - 0.5*lambda_star*(hB1_r- hB1_l)
    f_hB2 = 0.5 * (fleft_hB2 + fright_hB2) - 0.5*lambda_star*( hB2_r- hB2_l)

    flux[0]   = f_h * mesure
    flux[1]   = f_hu * mesure
    flux[2]   = f_hv * mesure
    flux[3]   = f_hB1 * mesure
    flux[4]   = f_hB2 * mesure
    flux[5]   = 0.
    flux[6]   = 0.




@inline
def srnh_scheme_MHD(hu_l:'float', hu_r:'float', hv_l:'float', hv_r:'float', h_l:'float', h_r:'float', hB1_l:'float', hB1_r:'float', hB2_l:'float',hB2_r:'float',
                hPSI_l:'float', hPSI_r:'float', hB1c_l:'float', hB1c_r:'float', hB2c_l:'float',hB2c_r:'float',
                Z_l:'float', Z_r:'float', normal:'float[:]', mesure:'float', grav:'float', flux:'float[:]', cpsi:'float'):
    
    from numpy import zeros, sqrt, fabs
    
    ninv =  zeros(2)
    w_dif =  zeros(6)

    ninv[0] = -1*normal[1]
    ninv[1] = normal[0]

    u_h = (hu_l / h_l * sqrt(h_l)
       + hu_r / h_r * sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))

    v_h = (hv_l / h_l * sqrt(h_l)
           + hv_r / h_r * sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))

    B1_h = (hB1_l / h_l * sqrt(h_l)
           + hB1_r / h_r * sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))

    B2_h = (hB2_l / h_l * sqrt(h_l)
           + hB2_r / h_r * sqrt(h_r)) /(sqrt(h_l) + sqrt(h_r))
        
    #uvh =  array([uh, vh])
    un_h = u_h*normal[0] + v_h*normal[1]
    un_h = un_h / mesure
    vn_h = u_h*ninv[0] + v_h*ninv[1]
    vn_h = vn_h / mesure

    B1n_h = B1_h*normal[0] + B2_h*normal[1]
    B1n_h = B1n_h / mesure
    B2n_h = B1_h*ninv[0] + B2_h*ninv[1]
    B2n_h = B2n_h / mesure

    hroe = (h_l+h_r)/2
    uroe = un_h
    vroe = vn_h
    B1roe = B1n_h
    B2roe = B2n_h
    
    uleft = hu_l*normal[0] + hv_l*normal[1]
    uleft = uleft / mesure
    vleft = hu_l*ninv[0] + hv_l*ninv[1]
    vleft = vleft / mesure

    B1left = hB1_l*normal[0] + hB2_l*normal[1]
    B1left = B1left / mesure
    B2left = hB1_l*ninv[0] + hB2_l*ninv[1]
    B2left = B2left / mesure

    uright = hu_r*normal[0] + hv_r*normal[1]
    uright = uright / mesure
    vright = hu_r*ninv[0] + hv_r*ninv[1]
    vright = vright / mesure

    B1right = hB1_r*normal[0] + hB2_r*normal[1]
    B1right = B1right / mesure
    B2right = hB1_r*ninv[0] + hB2_r*ninv[1]
    B2right = B2right / mesure


    B1i = (hB1c_l*normal[0] + hB2c_l*normal[1]) / mesure
    B1j = (hB1c_r*normal[0] + hB2c_r*normal[1])/ mesure     
                       
    absy = (B1i + B1j)/2 

    w_lrh = (h_l  + h_r)/2
    w_lrhu = (uleft + uright)/2
    w_lrhv = (vleft + vright)/2
    w_lrhB1 = (B1left + B1right)/2
    w_lrhB2 = (B2left + B2right)/2
    w_lrhPSI = (hPSI_l + hPSI_r)/2
    w_lrz = (Z_l + Z_r)/2
    
    w_dif[0] = h_r - h_l
    w_dif[1] = uright - uleft
    w_dif[2] = vright - vleft
    w_dif[3] = B1right  - B1left
    w_dif[4] = B2right  - B2left    
    w_dif[5] = Z_r - Z_l
    

    signA= zeros((6, 6))
        
    sound = sqrt(grav * hroe)
    
   # B1roe = absy/hroe
    
    w=sqrt(B1roe * B1roe + grav * hroe)

    
    lambda1 = uroe - w
    
    lambda2 = uroe - B1roe
    
    lambda3 = uroe + B1roe
    
    lambda4 = uroe + w
    

    epsilon = 1e-15

    if fabs(lambda1) < epsilon:
        s1 = 0.
        pi1 = 0.
    else:
        s1 = lambda1 /  fabs(lambda1)
        pi1 = s1/lambda1
        
    if  fabs(lambda2) < epsilon:
        s2 = 0.
        pi2 = 0.
    else:
        s2 = lambda2 /  fabs(lambda2)
        pi2 = 1./fabs(lambda2) 
    
    if   fabs(lambda3) < epsilon:
        s3 = 0.
        pi3 = 0.
    else:
        s3 = lambda3 /  fabs(lambda3)
        pi3 = 1./fabs(lambda3)
    
    if  fabs(lambda4) < epsilon:
        s4 = 0.
        pi4 = 0.
    else:
        s4 = lambda4 /  fabs(lambda4)
        pi4 = 1./fabs(lambda4)

    gamma1 = vroe + B2roe
    gamma2 = vroe - B2roe
         
    sigma1 = vroe*(s1*lambda4  - s4*lambda1)  - w*(s2*gamma1 + s3*gamma2)
    sigma2 = B2roe*(s1*lambda4 - s4*lambda1)  - w*(s2*gamma1 - s3*gamma2)
    
    if fabs(lambda2) < epsilon and fabs(lambda3) < epsilon  :

            mu1 = B1roe*vroe*pi1/w  - B1roe*vroe*pi4/w
            mu2 = B1roe*B2roe*pi1/w - B1roe*B2roe*pi4/w
            ann = 1         

    else :   
            mu1 = B1roe*vroe*pi1/w  - B1roe*vroe*pi4/w  - 0.5*(gamma1*pi2 - gamma2*pi3)  
            mu2 = B1roe*B2roe*pi1/w - B1roe*B2roe*pi4/w - 0.5*(gamma1*pi2 +  gamma2*pi3) 
            ann = 1

        
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
    w_lrhP  = w_lrhPSI  - u_hP
    w_lrhB1 = w_lrhB1 - u_hB1 #absy#w_lrhB1 - u_hB1
    w_lrhB2 = w_lrhB2 - u_hB2
    w_lrz   = w_lrz   - u_z
    


    w_hP   = hPSI_r - hPSI_l
    mw_hB1 = w_lrhB1/2 - (1/(2*cpsi))*w_hP
    mhP   = (hPSI_r - hPSI_l)/2 - cpsi*w_dif[3]/2
    

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
    
    norm = normal/mesure

    q_s = normal[0] * unew  + normal[1] * vnew
    p_s = normal[0] * B1new + normal[1] * B2new
    
    Flux_B1psi = mhP * norm[0]*mesure
    Flux_B2psi = mhP * norm[1]*mesure
    Flux_hPpsi = cpsi*cpsi*mw_hB1*mesure
    


    flux[0]   = q_s
    flux[1]   = q_s * w_lrhu/w_lrh + 0.5 * grav * w_lrh * w_lrh * normal[0] - p_s*w_lrhB1/w_lrh
    flux[2]   = q_s * w_lrhv/w_lrh + 0.5 * grav * w_lrh * w_lrh * normal[1] - p_s*w_lrhB2/w_lrh
    flux[3]   = (w_lrhv*w_lrhB1/w_lrh - w_lrhu*w_lrhB2/w_lrh ) * normal[1]  + Flux_B1psi 
    flux[4]   = (w_lrhu*w_lrhB2/w_lrh - w_lrhv*w_lrhB1/w_lrh ) * normal[0]  + Flux_B2psi  
    flux[5]   = Flux_hPpsi
    flux[6]   = 0.






def explicitscheme_convective_SWMHD(rez_h:'float[:]', rez_hu:'float[:]', rez_hv:'float[:]', rez_hB1:'float[:]', rez_hB2:'float[:]', rez_PSI:'float[:]', rez_Z:'float[:]', 
                                 h_c:'float[:]', hu_c:'float[:]', hv_c:'float[:]', hB1_c:'float[:]', hB2_c:'float[:]', hPSIc:'float[:]', Z_c:'float[:]', 
                                 h_ghost:'float[:]', hu_ghost:'float[:]', hv_ghost:'float[:]', hB1_ghost:'float[:]', hB2_ghost:'float[:]', hPSIghost:'float[:]', Z_ghost:'float[:]',
                                 h_halo:'float[:]', hu_halo:'float[:]', hv_halo:'float[:]', hB1_halo:'float[:]', hB2_halo:'float[:]', hPSIhalo:'float[:]', Z_halo:'float[:]',
                                 h_x:'float[:]', h_y:'float[:]', hx_halo:'float[:]', hy_halo:'float[:]',
                                 psi:'float[:]', psi_halo:'float[:]', 
                                 centerc:'float[:,:]', centerf:'float[:,:]', centerh:'float[:,:]', centerg:'float[:,:]',
                                 cellidf:'int[:,:]', mesuref:'float[:]', normalf:'float[:,:]', halofid:'int[:]',
                                 innerfaces:'int[:]', halofaces:'int[:]', boundaryfaces:'int[:]', periodicboundaryfaces:'int[:]', shift:'float[:,:]',
                                 order:'int', cpsi:'float', hB1_cst:'float[:]', hB2_cst:'float[:]'):


    rez_h[:] = 0.; rez_hu[:] = 0.; rez_hv[:] = 0.; rez_hB1[:] = 0.; rez_hB2[:] = 0.; rez_PSI[:] = 0.; rez_Z[:] = 0.

    
    from numpy import zeros
    
    grav = 1.0
    flux = zeros(7)
#    r_l = zeros(2)
#    r_r = zeros(2)
#    


    
    
    for i in innerfaces:
        
        h_l   = h_c[cellidf[i][0]]
        hu_l  = hu_c[cellidf[i][0]]
        hv_l  = hv_c[cellidf[i][0]]
        hB1_l = hB1_c[cellidf[i][0]]
        hB2_l = hB2_c[cellidf[i][0]]
        hPSI_l = hPSIc[cellidf[i][0]]
        Z_l   = Z_c[cellidf[i][0]]
        
        hB1c_l = hB1_cst[cellidf[i][0]]
        hB2c_l = hB2_cst[cellidf[i][0]]


        normal = normalf[i]
        mesure = mesuref[i]

        h_r   = h_c[cellidf[i][1]]
        hu_r  = hu_c[cellidf[i][1]]
        hv_r  = hv_c[cellidf[i][1]]
        hB1_r = hB1_c[cellidf[i][1]]
        hB2_r = hB2_c[cellidf[i][1]]
        hPSI_r = hPSIc[cellidf[i][1]]
        Z_r   = Z_c[cellidf[i][1]]
        hB1c_r = hB1_cst[cellidf[i][1]]
        hB2c_r = hB2_cst[cellidf[i][1]]            
        
#        center_left = centerc[cellidf[i][0]]
#        center_right = centerc[cellidf[i][1]]
#
#        h_x_left = h_x[cellidf[i][0]];   h_x_right = h_x[cellidf[i][1]]
#        h_y_left = h_y[cellidf[i][0]];   h_y_right = h_y[cellidf[i][1]]
#      
#        psi_left = psi[cellidf[i][0]]; psi_right = psi[cellidf[i][1]]
#        
#        r_l[0] = centerf[i][0] - center_left[0]; r_r[0] = centerf[i][0] - center_right[0]; 
#        r_l[1] = centerf[i][1] - center_left[1]; r_r[1] = centerf[i][1] - center_right[1]; 
        
        #h_l  = h_l  + (order - 1) * psi_left  * (h_x_left * r_l[0]  + h_y_left * r_l[1] )
        #h_r  = h_r  + (order - 1) * psi_right * (h_x_right* r_r[0]  + h_y_right* r_r[1] )
        
  
        srnh_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, hPSI_l, hPSI_r, hB1c_l, hB1c_r, hB2c_l, hB2c_r, Z_l, Z_r, normal, mesure, grav, flux,cpsi)
        #LF_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, normal, mesure, grav, flux)
            
        rez_h[cellidf[i][0]]  -= flux[0]
        rez_hu[cellidf[i][0]] -= flux[1]
        rez_hv[cellidf[i][0]] -= flux[2]
        rez_hB1[cellidf[i][0]] -= flux[3]
        rez_hB2[cellidf[i][0]] -= flux[4]
        rez_PSI[cellidf[i][0]] -= flux[5]                
        rez_Z[cellidf[i][0]]  -= flux[6]
          
        rez_h[cellidf[i][1]]  += flux[0]
        rez_hu[cellidf[i][1]] += flux[1]
        rez_hv[cellidf[i][1]] += flux[2]
        rez_hB1[cellidf[i][1]] += flux[3]
        rez_hB2[cellidf[i][1]] += flux[4]
        rez_PSI[cellidf[i][1]] += flux[5]
        rez_Z[cellidf[i][1]]  += flux[6]

    for i in halofaces:
        
        h_l   = h_c[cellidf[i][0]]
        hu_l  = hu_c[cellidf[i][0]]
        hv_l  = hv_c[cellidf[i][0]]
        hB1_l = hB1_c[cellidf[i][0]]
        hB2_l = hB2_c[cellidf[i][0]]
        hPSI_l = hPSIc[cellidf[i][0]]
        Z_l   = Z_c[cellidf[i][0]]
        hB1c_l = hB1_cst[cellidf[i][0]]
        hB2c_l = hB2_cst[cellidf[i][0]]



        normal = normalf[i]
        mesure = mesuref[i]
        h_r  = h_halo[halofid[i]]
        hu_r = hu_halo[halofid[i]]
        hv_r = hv_halo[halofid[i]]
        hB1_r = hB1_halo[halofid[i]]
        hB2_r = hB2_halo[halofid[i]]
        hPSI_r = hPSIhalo[halofid[i]]
        Z_r  = Z_halo[halofid[i]]
 
        hB1c_r = hB1_cst[cellidf[i][1]]
        hB2c_r = hB2_cst[cellidf[i][1]]            
         
#        
#        center_left = centerc[cellidf[i][0]]
#        center_right = centerh[halofid[i]]
#
#        h_x_left = h_x[cellidf[i][0]];   h_x_right = hx_halo[halofid[i]]
#        h_y_left = h_y[cellidf[i][0]];   h_y_right = hy_halo[halofid[i]]
#
#        psi_left = psi[cellidf[i][0]]; psi_right = psi_halo[halofid[i]]
#        
#        r_l[0] = centerf[i][0] - center_left[0]; r_r[0] = centerf[i][0] - center_right[0]; 
#        r_l[1] = centerf[i][1] - center_left[1]; r_r[1] = centerf[i][1] - center_right[1]; 
        
        #h_l  = h_l  + (order - 1) * psi_left  * (h_x_left   * r_l[0] + h_y_left   * r_l[1] )
        #h_r  = h_r  + (order - 1) * psi_right * (h_x_right  * r_r[0] + h_y_right  * r_r[1] )
       
        
        srnh_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, hPSI_l, hPSI_r, hB1c_l, hB1c_r, hB2c_l, hB2c_r, Z_l, Z_r, normal, mesure, grav, flux, cpsi)

        #LF_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, normal, mesure, grav, flux)


        rez_h[cellidf[i][0]]   -= flux[0]
        rez_hu[cellidf[i][0]]  -= flux[1]
        rez_hv[cellidf[i][0]]  -= flux[2]
        rez_hB1[cellidf[i][0]] -= flux[3]
        rez_hB2[cellidf[i][0]] -= flux[4]
        rez_PSI[cellidf[i][0]] -= flux[5]                
        rez_Z[cellidf[i][0]]   -= flux[6]            


    for i in periodicboundaryfaces:
        

        h_l   = h_c[cellidf[i][0]]
        hu_l  = hu_c[cellidf[i][0]]
        hv_l  = hv_c[cellidf[i][0]]
        hB1_l = hB1_c[cellidf[i][0]]
        hB2_l = hB2_c[cellidf[i][0]]
        hPSI_l = hPSIc[cellidf[i][0]]
        Z_l   = Z_c[cellidf[i][0]]
        
        hB1c_l = hB1_cst[cellidf[i][0]]
        hB2c_l = hB2_cst[cellidf[i][0]]


        normal = normalf[i]
        mesure = mesuref[i]

        h_r   = h_c[cellidf[i][1]]
        hu_r  = hu_c[cellidf[i][1]]
        hv_r  = hv_c[cellidf[i][1]]
        hB1_r = hB1_c[cellidf[i][1]]
        hB2_r = hB2_c[cellidf[i][1]]
        hPSI_r = hPSIc[cellidf[i][1]]
        Z_r   = Z_c[cellidf[i][1]]
        hB1c_r = hB1_cst[cellidf[i][1]]
        hB2c_r = hB2_cst[cellidf[i][1]]
        
        srnh_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, hPSI_l, hPSI_r, hB1c_l, hB1c_r, hB2c_l, hB2c_r, Z_l, Z_r, normal, mesure, grav, flux, cpsi)
        
        rez_h[cellidf[i][0]]   -= flux[0]
        rez_hu[cellidf[i][0]]  -= flux[1]
        rez_hv[cellidf[i][0]]  -= flux[2]
        rez_hB1[cellidf[i][0]] -= flux[3]
        rez_hB2[cellidf[i][0]] -= flux[4]
        rez_PSI[cellidf[i][0]] -= flux[5]                
        rez_Z[cellidf[i][0]]   -= flux[6]

        
    
    for i in boundaryfaces:
       
        h_l   = h_c[cellidf[i][0]]
        hu_l  = hu_c[cellidf[i][0]]
        hv_l  = hv_c[cellidf[i][0]]
        hB1_l = hB1_c[cellidf[i][0]]
        hB2_l = hB2_c[cellidf[i][0]]
        hPSI_l = hPSIc[cellidf[i][0]]
        Z_l   = Z_c[cellidf[i][0]]
        hB1c_l = hB1_cst[cellidf[i][0]]
        hB2c_l = hB2_cst[cellidf[i][0]]

   
        normal = normalf[i]
        mesure = mesuref[i]
        h_r  = h_ghost[i]
        hu_r = hu_ghost[i]
        hv_r = hv_ghost[i]
        hB1_r = hB1_ghost[i]
        hB2_r = hB2_ghost[i]
        hPSI_r = hPSIghost[i]
        Z_r  = Z_ghost[i]

        hB1c_r = hB1_cst[cellidf[i][1]]
        hB2c_r = hB2_cst[cellidf[i][1]]            
 
        
#        center_left = centerc[cellidf[i][0]]
#        center_right = centerg[i]
#
#        h_x_left = h_x[cellidf[i][0]];   h_y_left = h_y[cellidf[i][0]]; 
#
#       
#        psi_left = psi[cellidf[i][0]]; 
#        
#        r_l[0] = centerf[i][0] - center_left[0]; r_r[0] = centerf[i][0] - center_right[0]; 
#        r_l[1] = centerf[i][1] - center_left[1]; r_r[1] = centerf[i][1] - center_right[1]; 
        
        #h_l  = h_l  + (order - 1) * psi_left  * (h_x_left * r_l[0] + h_y_left * r_l[1] )
        #h_r  = h_r 


        
        srnh_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, hPSI_l, hPSI_r, hB1c_l, hB1c_r, hB2c_l, hB2c_r, Z_l, Z_r, normal, mesure, grav, flux, cpsi)
        #LF_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, normal, mesure, grav, flux)
            
        rez_h[cellidf[i][0]]   -= flux[0]
        rez_hu[cellidf[i][0]]  -= flux[1]
        rez_hv[cellidf[i][0]]  -= flux[2]
        rez_hB1[cellidf[i][0]] -= flux[3]
        rez_hB2[cellidf[i][0]] -= flux[4]
        rez_PSI[cellidf[i][0]] -= flux[5]                
        rez_Z[cellidf[i][0]]   -= flux[6]            
