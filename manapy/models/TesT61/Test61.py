#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 09:16:33 2022

@author: cisse
"""

from manapy.ddm1 import Domain, readmesh

from manapy.models.SWMHDModel import ShallowWaterMHDModel, initialisation_SWMHD, Total_Energy

from manapy.ast import Variable

import numpy as np
import timeit
import os
from mpi4py import MPI


          
    
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

#File name

filename = "EC3.msh"

filename = os.path.join(MESH_DIR, filename)
dim = 2
readmesh(filename, dim=dim, periodic=[0,0,0])

#Create the informations about cells, faces and nodes
domain = Domain(dim=dim)

faces = domain.faces
cells = domain.cells
halos = domain.halos
nodes = domain.nodes

nbcells = domain.nbcells

if RANK == 0: print("Start Computation")

#TODO tfinal
choix = 1
grav = 1.0
cfl = 0.8
time = 0
tfinal = 0.2
order = 1
saving_at_node = 1

# wave vectors
k1 = 2.0 ; k2 = 0.0 ; eps = .1 ; tol = .01 #0.1


parameters = {}
parameters["Manning"] = 0.
parameters["CoriolisForce"] = 0.
parameters["vepsilon"] = .0



SWMHDModel = ShallowWaterMHDModel(domain, terms=['source', 'stabilizater'], parameters=parameters, order=order)

h   = SWMHDModel.h
hu  = SWMHDModel.hu
hv  = SWMHDModel.hv
hB1 = SWMHDModel.hB1
hB2 = SWMHDModel.hB2
PSI = SWMHDModel.PSI
Z   = SWMHDModel.Z

g = 1.0  
            



boundaries = {"in" : "nonslip",
              "out" : "nonslip",
              "upper":"nonslip",
              "bottom":"nonslip"
              }
values = {}

#SWMHDModel.h   = Variable(domain=domain, BC=boundaries, terms=['convective', 'source', 'coriolis'], values=values)

SWMHDModel.hu  = Variable(domain=domain, BC=boundaries, terms=['convective', 'source',  'coriolis'], values=values)
SWMHDModel.hv  = Variable(domain=domain, BC=boundaries, terms=['convective', 'source',  'coriolis'], values=values)
SWMHDModel.hB1 = Variable(domain=domain, BC=boundaries, terms=['convective', 'source',  'coriolis'], values=values)
SWMHDModel.hB2 = Variable(domain=domain, BC=boundaries, terms=['convective', 'source',  'coriolis'], values=values)
SWMHDModel.PSI = Variable(domain=domain, BC=boundaries, terms=['convective', 'source',  'coriolis'], values=values)



#SWMHDModel.PSI = Variable(domain=domain, BC=boundaries, terms=['coriolis'], values=values)


####Initialisation

# import sys;
initialisation_SWMHD(h.cell, hu.cell, hv.cell, hB1.cell, hB2.cell, PSI.cell, Z.cell, cells.center, choix, k1, k2, eps, tol)

et, ec, em, ep = Total_Energy(h.cell, hu.cell, hv.cell, hB1.cell, hB2.cell, Z.cell, grav, cells._volume)
Et = []
Ec = []
Em = []
Ep = []
TimeS = []

TimeS.append(0.0)
Et.append(et)
Ec.append(ec)
Em.append(em)
Ep.append(ep)


#update time step
SWMHDModel.update_time_step(cfl=cfl)

#saving 50 vtk file
tot = int((tfinal-time)/SWMHDModel.time_step/50)
miter = 0
niter = 1

if RANK == 0:
    print("Start loop")
    
start = timeit.default_timer()

#loop over time
while time < tfinal:
    
    time = time + SWMHDModel.time_step
    
    #update the ghost values for the boundary conditions
    SWMHDModel.update_values()
    
    # update the global value    
    SWMHDModel.update_cpsi_global(cfl=cfl, comm=COMM)
        
    #update solution   
    h_n = h.cell
    #h_n.interpolate_celltonode()
#    hu_n= hu
#    hv_n= hv
#    hB1_n= hB1
#    hB2_n= hB2
    
    SWMHDModel.update_solution()


    hexact     = np.zeros(nbcells) 
    uexact     = np.zeros(nbcells) 
    vexact     = np.zeros(nbcells) 
    B1exact    = np.zeros(nbcells) 
    B2exact    = np.zeros(nbcells) 
    div        = np.zeros(nbcells)
    Gd         = np.zeros(nbcells) 

    

    #update time step
    
    SWMHDModel.update_time_step(cfl=cfl)
    
    h.interpolate_celltonode()
    hu.interpolate_celltonode()
    hv.interpolate_celltonode()
    hB1.interpolate_celltonode()
    hB2.interpolate_celltonode()
    Z.interpolate_celltonode()
    
    hB1.compute_cell_gradient()
    hB2.compute_cell_gradient()
    
    FFG = hB1.gradcellx + hB2.gradcelly


    
    



    
    f1 = lambda x, y, z : 100
    f2 = lambda x, y, z : 2 + np.sin(2*np.pi*(x - y))
    f3 = lambda x, y, z : 5 + 2*np.sin(2*np.pi*(x - y))
       
    for i in range(nbcells):
        
        xcent = cells.center[i][0]
        ycent = cells.center[i][1]
        #he[i] = f1(xcent, ycent, 0.)
        #ue[i] = f2(xcent, ycent, 0.)
        #Be[i] = f3(xcent, ycent, 0.)
        div[i] = FFG[i]*cells._volume[i]
        
        
    #save vtk files for the solution
    if niter== 1 or niter%tot == 0:
        if saving_at_node:
            #print("Error L1 in h =", hn.norml2(h_n, order=1))
            et, ec, em, ep = Total_Energy(h.node, hu.node, hv.node, hB1.node, hB2.node, Z.node, grav, cells._volume)
            TimeS.append(time)
            Et.append(et)
            Ec.append(ec)
            Em.append(em)
            Ep.append(ep)
            
            #print("l1 norm for h is ", errorh) 
            #print("l1 norm for u is ", erroru)
            #print("l1 norm for B is ", errorB)   
            
            domain.save_on_node(SWMHDModel.time_step, time, niter, miter, valueh=h.node, valuehu=hu.node, valuehv=hv.node,
                                valuehB1=hB1.node, valuehB2=hB2.node, valueZ=Z.node)



        else:
            domain.save_on_cell(SWMHDModel.time_step, time, niter, miter, value=h.cell, valuex=hexact)
        miter += 1

    niter += 1

if saving_at_node:
    #print("Error L1 in h =", hn.norml2(h_n, order=1))
    et, ec, em, ep = Total_Energy(h.node, hu.node, hv.node, hB1.node, hB2.node, Z.node, grav, cells._volume)
    TimeS.append(time)
    Et.append(et)
    Ec.append(ec)
    Em.append(em)
    Ep.append(ep)
               
    domain.save_on_node(SWMHDModel.time_step, time, niter, miter, valueh=h.node, valuehu=hu.node, valuehv=hv.node,
                        valuehB1=hB1.node, valuehB2=hB2.node, valueZ=Z.node)
else:

    domain.save_on_cell(SWMHDModel.time_step, time, niter, miter, value=h.cell, valuex=hexact)

#tri_plot(h.node, hu.node, nodes._cellid, nodes._vertex, time)
with open("Temps.txt", "w") as f :
      np.savetxt(f, TimeS, fmt='%.15f')
with open("ET.txt", "w") as f :
      np.savetxt(f, Et, fmt='%.15f')
with open("EC.txt", "w") as f :
      np.savetxt(f, Ec, fmt='%.15f')
with open("EM.txt", "w") as f :
      np.savetxt(f, Em, fmt='%.15f')            
with open("Ep.txt", "w") as f :
      np.savetxt(f, Ep, fmt='%.15f')


stop = timeit.default_timer()

if RANK == 0: print(stop - start)
