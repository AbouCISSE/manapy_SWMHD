#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 20:51:37 2020

@author: kissami
"""
import timeit
import os
import numpy as np
from mpi4py import MPI
import manapy.ddm as ddm 
import manapy.ddp as ddp
from manapy.ddp                import readmesh, generate_structure
from manapy.comms              import (create_mpi_graph, prepare_comm, define_halosend, 
                                       all_to_all, update_haloghost_info_2d)


COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

# ... get the mesh directory
try:
    MESH_DIR = os.environ['MESH_DIR']

except:
    BASE_DIR = os.path.dirname(os.path.realpath(__file__))
    BASE_DIR = os.path.join(BASE_DIR , '..', '..', '..')
    MESH_DIR = os.path.join(BASE_DIR, 'mesh')
#TODO mesh
dim = 2
#File name
if dim == 2:
 
    filename = os.path.join(MESH_DIR, "EC3.msh")
 
    typeofCells = 'triangle'
elif dim == 3:
    filename = os.path.join(MESH_DIR, "cube.msh")
    typeofCells = 'tetra'
print(RANK, SIZE)
if RANK == 0:
    #Read gmsh file
    readmesh(filename, dim=dim, size=SIZE, periodic=[0,1])
    #removing existing vtk files
    mypath = "results"
    if not os.path.exists(mypath):
        os.mkdir(mypath)
    for root, dirs, files in os.walk(mypath):
        for file in files:
            os.remove(os.path.join(root, file))

COMM.Barrier()
#Create the informations about cells, faces and nodes
grid = generate_structure(dim=dim, size=SIZE)


faces = grid["faces"]
cells = grid["cells"]
halos = grid["halos"]
nodes = grid["nodes"]

nbelements = len(cells.center)
nbfaces = len(faces.name)
nbnodes = len(nodes.vertex)

mystruct = np.dtype([('h',  np.double),
                     ('hu', np.double),
                     ('hv', np.double),
                     ('hP', np.double),
                     ('hB1', np.double),
                     ('hB2', np.double),
                     ('Ch', np.double),
                     ('Z',  np.double),])
variables = mystruct.names 

w_c     = np.zeros(nbelements, dtype=mystruct, order='C')
w_p     = np.zeros(nbfaces, dtype=mystruct, order='C')
w_x     = np.zeros(nbelements, dtype=mystruct, order='C')
w_y     = np.zeros(nbelements, dtype=mystruct, order='C')
psi     = np.zeros(nbelements, dtype=np.double)
src     = np.zeros(nbelements, dtype=mystruct, order='C')
rez     = np.zeros(nbelements, dtype=mystruct, order='C')
cor     = np.zeros(nbelements, dtype=mystruct, order='C')
w_ghost = np.zeros(nbfaces, dtype=mystruct, order='C')
wx_face = np.zeros(nbfaces, dtype=mystruct, order='C')
wy_face = np.zeros(nbfaces, dtype=mystruct, order='C')
wexact  = np.zeros(nbelements, dtype=mystruct, order='C')
w_node = np.zeros(nbnodes, dtype=mystruct, order='C')
wnode_exact = np.zeros(nbnodes, dtype=mystruct, order='C')

sizehaloghost = update_haloghost_info_2d(nodes, cells, halos)
#ghost of halo cells
w_haloghost = np.zeros(sizehaloghost,dtype=mystruct, order='C')

#compute the arrays needed for the mpi communication
scount, sdepl, rcount, rdepl, taille, indsend = ddm.prepare_comm(cells, halos)

#compute the interpolation variables
R_x, R_y, lambda_x,lambda_y, number =  ddm.variables(cells.center, nodes.cellid, nodes.halonid, nodes.vertex,
                                                     nodes.name, halos.centvol, nodes.ghostcenter, 
                                                     nodes.haloghostcenter,cells.shift,nodes.periodicid)
#R_x, R_y, lambda_x,lambda_y, number =  ddp.variablesNP(cells.center, nodes.cellid, nodes.halonid,
#                                                     nodes.vertex, nodes.name, nodes.ghostcenter,
#                                                     halos.centvol)
w_halosend   = np.zeros(len(halos.halosint), dtype=mystruct)
wx_halosend  = np.zeros(len(halos.halosint), dtype=mystruct)
wy_halosend  = np.zeros(len(halos.halosint), dtype=mystruct)
psi_halosend = np.zeros(len(halos.halosint), dtype=np.double)

w_halo   = np.zeros(taille, dtype=mystruct)
wx_halo  = np.zeros(taille, dtype=mystruct)
wy_halo  = np.zeros(taille, dtype=mystruct)
psi_halo = np.zeros(taille, dtype=np.double)
print("sizeshift =", len(cells.shift))
w_n = w_c
#=-=*=-=*=-=*=-=*=-=*=-=*=-=*=-=*=-=*=-=*Initialisation=-=*=-=*=-=*=-=*=-=*=-=*=-=*=-=*=-=*=-=*#

choix        = 2

w_c = ddm.initialisation(w_c, cells.center, choix)
w_w = ddm.initialisation(w_c, cells.center, choix)
    

    



#update the halo values
w_halosend = np.zeros(len(halos.halosint), dtype=mystruct)

#TODO tfinal
if RANK == 0: print("Start Computation")
cfl = 0.8
time = 0
tfinal = 0.4
iterfinal = 9000000
SCHEME =  "SRNH"  # "LF" or "GLM" or "SRNH-GLM" or "SRNH"
order = 1
saving_at_node = 1

term_convec   = "on"
term_source   = "off"
term_coriolis = "off"
calcul_div    = "on"
globav        = "on"
MGLM          = "off"  # if SCHEME =  "SRNH" or "LF" # "on" or "off" or "crt"  
#MGLM          = "on"  # if SCHEME =  "SRNH-GLM" or "GLM" 
#MGLM          = "crt" # if SCHEME =  "SRNH-GLM" or "GLM" 
correction    = "off" 
grav = 1.00
cd   = 0.14
cr   = .01# 0.01
delta = 1
f_c = 0#0.5
TimeS = []
L2div = []
L1div = []
Linfdiv = []
ENEGT = []
TimeS.append(0.0)
L2div.append(0.0)
L1div.append(0.0)
Linfdiv.append(0.0)
div0, E_t = ddm.norml2(w_c, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ddm.Total_Energy(w_c, grav,correction))
ENEGT.append(E_t)   
div=0.0
print("******************delta=", delta,"*********************")
print("******************cd=", cd,"*********************")

###calculation of the time step
d_t = ddm.time_step(w_c, cfl, faces.normal, faces.mesure, cells.volume, cells.faceid, grav, term_convec)
dt_i = np.zeros(1)
COMM.Allreduce(d_t, dt_i, MPI.MIN)
d_t =  np.double(dt_i)

#saving 50 vtk file
tot = int(tfinal/d_t/50)
if tot==0:
    tot =1
miter = 0
niter = 1


if RANK == 0: 
    print("Start loop")
timec=0
w_ghostc = ddm.ghost_value(w_c, w_ghost, faces.cellid, faces.name, faces.normal, faces.mesure,
                              timec, cells.center, faces.center, cells.cellnid, w_x, w_y)
start = timeit.default_timer()
#print(cells.cellnid)
#loop over time
#c_h = ddm.global_cpsi(w_c, cfl, faces.normal, faces.mesure, cells.volume, cells.faceid, grav, globav)
#c_h = c_h*0.01

while time < tfinal and niter < iterfinal:

    time = time + d_t
    
    c_h = ddm.global_cpsi(w_c, cfl, faces.normal, faces.mesure, cells.volume, cells.faceid, grav, globav)


    #update the ghost values for the boundary conditions
    w_ghost = ddm.ghost_value(w_c, w_ghost, faces.cellid, faces.name, faces.normal, faces.mesure,
                              time, cells.center, faces.center, cells.cellnid, w_x, w_y)
    
   #update the halo values
    w_halosend = ddm.define_halosend(w_c, w_halosend, indsend)
    ddm.all_to_all_struct(w_halosend, taille, mystruct, variables,
                          scount, sdepl, rcount, rdepl, w_halo)
#    w_halosend = ddp.define_halosend(w_c, w_halosend, indsend, SIZE)
#    ddp.all_to_all(w_halosend, taille, scount, sdepl, rcount, rdepl, w_halo, comm_ptr)
    
    w_haloghost = ddm.update_haloghost_value_MHD(w_halo, w_haloghost, nodes.haloghostcenter, 
                                                 nodes.vertex, nodes.name, time)

    

    if calcul_div == "on" :
        #compute derivative
        w_x, w_y = ddp.cell_gradient2(w_c, w_ghost, w_halo, w_haloghost, cells.center, nodes.vertex, cells.nodeid, 	
                                      cells.cellnid, cells.halonid, nodes.haloghostcenter, nodes.name, nodes.ghostcenter, 
                                      halos.centvol, w_x, w_y, cells.shift, cells.periodicnid)
                                                         
        
    
#        rez_div = ddm.explicitscheme_div(w_c, w_p, w_x, w_y, psi, w_ghost, w_halo, wx_halo, wy_halo, psi_halo,
#                           faces.cellid, cells.faceid, cells.nodeid, cells.center, cells.cellfid, halos.centvol, 
#                           faces.mesure, faces.center, faces.normal, faces.halofid, faces.name,
#                           faces.ghostcenter, nodes.cellid, cells.shift, mystruct, order, SCHEME, grav)
        
        
        DIV = np.zeros(len(w_c))#ddm.calculate_div(w_c, rez_div, cells.volume)
        
        #update the halo  derivatives values
        #wx_halosend = ddm.define_halosend(w_x, wx_halosend, indsend)
        #wy_halosend = ddm.define_halosend(w_y, wy_halosend, indsend)

        #ddm.all_to_all_struct(wx_halosend, taille, mystruct, variables,
        #                      scount, sdepl, rcount, rdepl, wx_halo)
        #ddm.all_to_all_struct(wy_halosend, taille, mystruct, variables,
        #                      scount, sdepl, rcount, rdepl, wy_halo)
    

   
    if order == 2 :
        #compute derivative
        w_x, w_y = ddp.cell_gradient2(w_c, w_ghost, w_halo, w_haloghost, cells.center, nodes.vertex, cells.nodeid, 	
                                      cells.cellnid, cells.halonid, nodes.haloghostcenter, nodes.name, nodes.ghostcenter, 
                                      halos.centvol, w_x, w_y, cells.shift, cells.periodicnid)
        
        ddm.barthlimiter(w_c, w_ghost, w_halo, w_x, w_y, psi, faces.cellid, cells.faceid, faces.name, 
                         faces.halofid, cells.center, faces.center)
        
        #update the halo  derivatives values
        wx_halosend = ddm.define_halosend(w_x, wx_halosend, indsend)
        wy_halosend = ddm.define_halosend(w_y, wy_halosend, indsend)
        psi_halosend = ddm.define_halosend(psi, psi_halosend, indsend)

        ddm.all_to_all_struct(wx_halosend, taille, mystruct, variables,
                              scount, sdepl, rcount, rdepl, wx_halo)
        ddm.all_to_all_struct(wy_halosend, taille, mystruct, variables,
                              scount, sdepl, rcount, rdepl, wy_halo)
        
        ddm.all_to_all(psi_halosend, taille, scount, sdepl, rcount, rdepl, 
                       psi_halo)
    
   
    w_node = np.zeros(nbnodes, dtype=mystruct, order='C')
    wnode_exact = np.zeros(nbnodes, dtype=mystruct, order='C')
    
    wexact = ddm.exact_solution(wexact, cells.center, time, grav, choix)
    w_node, wnode_exact = ddm.centertovertex(w_c, wexact, w_ghost, w_halo, w_haloghost, nodes.haloghostcenter, 
                                             cells.center, halos.centvol, nodes.cellid, nodes.halonid, nodes.vertex,  
                                             nodes.name, nodes.ghostcenter, w_node, wnode_exact, R_x, R_y, lambda_x, 
                                             lambda_y, number, cells.shift, nodes.periodicid)

       
    
    #ddm.plot_tricontour( w_node["h"], nodes.cellid, cells.center, nodes.vertex)
    
    if (term_convec == "on"):
        #update the rezidus using explicit scheme
        rez = ddm.explicitscheme_convective(w_c, w_p, w_x, w_y, psi, w_ghost, w_halo, wx_halo, wy_halo, psi_halo,
                                           faces.cellid, cells.faceid, cells.nodeid, cells.center, 
                                           cells.cellfid, halos.centvol, faces.mesure,
                                           faces.center, faces.normal, faces.halofid, faces.name, 
                                           faces.ghostcenter, nodes.cellid, cells.shift, mystruct, order, SCHEME, grav, 
                                           c_h,delta, w_w, w_ghostc)

    


    if (term_source == "on"):
         #update the source using explicit scheme
         ddm.term_source_srnh_SW(src, w_c, w_ghost, w_halo, w_x, w_y, psi, wx_halo, wy_halo, psi_halo, 
                             cells.nodeid, cells.faceid, cells.cellfid, faces.cellid, cells.center,
                             cells.nf, faces.name, faces.center, halos.centvol,
                             nodes.vertex, faces.halofid, grav, order)



#        src = ddm.term_source_srnh(w_c, w_p, w_ghost, w_halo, w_x, w_y, psi, wx_halo, wy_halo, psi_halo,
#                                  cells.nodeid, cells.faceid, cells.cellfid, cells.center, cells.volume, 
#                                  cells.nf, faces.cellid, faces.nodeid, faces.normal, faces.center, 
#                                  faces.name, faces.ghostcenter, nodes.vertex, faces.halofid, halos.centvol,
#                                  mystruct, order, grav, src)
 
                                  
    w_n = ddm.update(w_c, w_n, d_t, rez, src, cor, cells.volume, cd, cr, c_h, MGLM)
    w_c = w_n
    #w_w = w_c


    #save vtk files for the solution
    if (niter == 1 or niter%tot == 0):
        if saving_at_node:

            ET = ddm.Total_Energy(w_node, grav, correction)
            ddm.save_paraview_results(w_node, wnode_exact, niter, miter, time, d_t, RANK, SIZE, cells.volume,
                                      cells.nodeid, nodes.vertex, w_x, w_y, tfinal, ET)
            
            div, ENT = ddm.norml2(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
            div1 = ddm.norml1(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
            divinf = ddm.normlinf(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
            TimeS.append(time)
            L2div.append(div)
            L1div.append(div1)
            Linfdiv.append(divinf)
            ENEGT.append(ENT)

       
        else:
            ET = ddm.Total_Energy(w_c, grav, correction)
            ddm.save_paraview_results(w_c, wexact, niter, miter, time, d_t, RANK, SIZE, cells.volume,
                                  cells.nodeid, nodes.vertex, w_x, w_y, tfinal, ET)
            div, ENT = ddm.norml2(w_c, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
            div1     = ddm.norml1(w_c, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
            divinf   = ddm.normlinf(w_c, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)


            TimeS.append(time)
            L2div.append(div)
            L1div.append(div1)
            Linfdiv.append(divinf)
            ENEGT.append(ENT)
        miter += 1

    
    niter += 1

    dt_c = np.zeros(nbelements)
    d_t =ddm.time_step(w_c, cfl, faces.normal, faces.mesure, cells.volume, cells.faceid, grav, term_convec)
    dt_i = np.zeros(1)
    COMM.Allreduce(d_t, dt_i, MPI.MIN)
    d_t = np.double(dt_i)
    tot = int(tfinal/d_t/50)
    if tot==0:
        tot =1

#ddm.plot_tricontour( w_node["h"], nodes.cellid, cells.center, nodes.vertex)
if saving_at_node:
    ET = ddm.Total_Energy(w_node, grav, correction)
    ddm.save_paraview_results(w_node, wnode_exact, niter, miter, time, d_t, RANK, SIZE, cells.volume,
                              cells.nodeid, nodes.vertex, w_x, w_y, tfinal, ET)
    div, ENT = ddm.norml2(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
    div1 = ddm.norml1(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
    divinf = ddm.normlinf(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)

    TimeS.append(time)
    L1div.append(div1)
    Linfdiv.append(divinf)
    L2div.append(div)       
    ENEGT.append(ENT)
else:

    ET = ddm.Total_Energy(w_c, grav, correction)
    ddm.save_paraview_results(w_c, wexact, niter, miter, time, d_t, RANK, SIZE, cells.volume,
                              cells.nodeid, nodes.vertex, w_x, w_y, tfinal, ET)
    div, ENT = ddm.norml2(w_c, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
    div1 = ddm.norml1(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)
    divinf = ddm.normlinf(w_node, cells.volume, cells.nodeid, nodes.vertex, w_x, w_y, ET)

    TimeS.append(time)
    L2div.append(div) 
    L1div.append(div1)
    Linfdiv.append(divinf)
    ENEGT.append(ENT)


stop = timeit.default_timer()

if RANK == 0: print(stop - start)

with open("Temps.txt", "w") as f :
      np.savetxt(f, TimeS, fmt='%.4f')
#with open("divhbL2.txt", "w") as f :
#      np.savetxt(f, L2div, fmt='%.4f')
#with open("divhbLinf.txt", "w") as f :
#      np.savetxt(f, Linfdiv, fmt='%.4f')
#with open("divhbL1.txt", "w") as f :
#      np.savetxt(f, L1div, fmt='%.4f')            
with open("ET.txt", "w") as f :
      np.savetxt(f, ENEGT, fmt='%.4f')

            


