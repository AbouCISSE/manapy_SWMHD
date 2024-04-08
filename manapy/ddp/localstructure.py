#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:08:27 2021

@author: kissami
"""

#from numba.typed import Dict
import timeit
from numpy import (zeros, ones, asarray,  double, int64, unique, where, 
                   array, sort, dot, cross, float32, append)
from collections import OrderedDict

from mpi4py import MPI

from manapy.ddp.module_ddp import (create_info_2dfaces, create_info_3dfaces,
                                   Compute_2dcentervolumeOfCell, Compute_3dcentervolumeOfCell,
                                   create_cellsOfFace, create_2dfaces, create_cell_faceid,
                                   create_3dfaces, create_NormalFacesOfCell)

from manapy.ddp.utils import (create_2doppNodeOfFaces, create_3doppNodeOfFaces,
                              create_NeighborCellByFace, create_node_cellid, 
                              oriente_3dfacenodeid)

COMM = MPI.COMM_WORLD
#SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


__all__ = ["generate_structure"]

class Cells:
    nodeid = []
    faceid = []
    center = []
    volume = []
    cellfid = []
    cellnid = []
    halonid = []
    nf = [] 
    globtoloc = OrderedDict()
    loctoglob = []
    tc = []
    shift = []
    periodicnid = []
    periodicfid = []
    
    
    def __init__(self, nodeid, faceid, center, volume, cellfid, cellnid, nf, globtoloc, 
                 loctoglob, halonid, tc, shift, periodicfid, periodicnid):
        self.nodeid = nodeid    # instance variable unique to each instance
        self.faceid = faceid
        self.center = center
        self.volume = volume
        self.cellfid = cellfid
        self.cellnid = cellnid
        self.halonid = halonid
        self.nf = nf
        self.globtoloc = globtoloc
        self.loctoglob = loctoglob
        self.tc = tc
        self.periodicnid = periodicnid
        self.periodicfid = periodicfid
        self.shift = shift
        
class Nodes:
    vertex = []
    name = []
    cellid = []
    halonid = []
    ghostcenter = []
    haloghostcenter = []
    loctoglob = []
    nparts = []
    periodicid = []

    def __init__(self, vertex, name, cellid, ghostcenter, loctoglob, halonid, nparts, 
                 haloghostcenter, periodicid):
        self.vertex = vertex    #instance variable unique to each instance
        self.name = name
        self.cellid = cellid
        self.ghostcenter = ghostcenter
        self.haloghostcenter = haloghostcenter
        self.loctoglob = loctoglob
        self.halonid = halonid
        self.nparts = nparts
        self.periodicid = periodicid
        
class Faces:
    nodeid = []
    cellid = []
    name = []
    normal = []
    mesure = []
    bound = 0
    center = []
    ghostcenter = []
    oppnodeid = []
    halofid = []

    def __init__(self, nodeid, cellid, name, normal, mesure, center, bound, ghostcenter, 
                 oppnodeid, halofid):
        self.nodeid = nodeid    # instance variable unique to each instance
        self.cellid = cellid
        self.name = name
        self.normal = normal
        self.mesure = mesure
        self.bound = bound
        self.center = center
        self.ghostcenter = ghostcenter
        self.oppnodeid = oppnodeid
        self.halofid = halofid

class Halo:
    halosint = []
    halosext = []
    neigh = []
    centvol = []
    faces = OrderedDict()
    nodes = OrderedDict()

    def __init__(self, halosint, halosext, centvol, neigh, faces, nodes):
        self.halosint = halosint    # instance variable unique to each instance
        self.halosext = halosext
        self.neigh = neigh
        self.centvol = centvol
        self.faces = faces
        self.nodes = nodes


def func_unique(array):
        uniq, index = unique(array, return_index=True)
        return uniq[index.argsort()]

def CreateStructure(file, dim, size):
    
    #Lecture des cellules à partir du fichier mesh..txt
    if size > 1 and RANK == 0:
        for line in file:
            if line == "GtoL\n":
                continue
            if line == "endGtoL\n":
                continue
            if line == "elements\n":
                break
            Cells.tc.append([int64(x) for x in line.split()][0])
    #Lecture des coordonnées des noeuds à partir du fichier mesh..txt
    for line in file:
        #read elements
        if line == "elements\n":
            continue
        if line == "endelements\n":
            continue
        if line == "nodes\n":
            break
        Cells.nodeid.append([int64(x) for x in line.split()])
    #Lecture des coordonnées des noeuds à partir du fichier mesh..txt
    for line in file:
        #read Nodes
        if line == "nodes\n":
            continue
        if line == "endnodes\n":
            continue
        if line == "halosint\n":
            break
        Nodes.vertex.append([double(x) for x in line.split()])
        
    Cells.nodeid = asarray(Cells.nodeid, dtype=int64)
    Nodes.vertex = asarray(Nodes.vertex, dtype=double)

    nbelements = len(Cells.nodeid)
    nbnodes    = len(Nodes.vertex)
    
    if size == 1:
        Cells.tc = zeros(nbelements, dtype=int64)
        #tc to gather solution of linear system
        for i in range(nbelements):
            Cells.tc[i] = i
    
    Nodes.name = zeros(nbnodes, dtype=int64)
    for i in range(nbnodes):
        Nodes.name[i] = int64(Nodes.vertex[i][3])
    
    #Create center and volume for each cell (pycceliser)
    Cells.center = zeros((nbelements, 3), dtype=double)
    Cells.volume = zeros(nbelements, dtype=double)
    if dim == 2:
        Compute_2dcentervolumeOfCell(Cells.nodeid, Nodes.vertex, nbelements, Cells.center, Cells.volume)
    elif dim == 3:
        Compute_3dcentervolumeOfCell(Cells.nodeid, Nodes.vertex, nbelements, Cells.center, Cells.volume)
    #create cells over each node (still numba function)
    Nodes.cellid, Cells.cellnid =  create_node_cellid(Cells.nodeid, Nodes.vertex, nbelements, nbnodes, dim=dim)
    
    #creating faces (pycceliser)
    faces = zeros(((dim+1)*nbelements, dim), dtype=int64)
    cellf = zeros((nbelements, dim+1), dtype=int64)
    
    if dim == 2:
        create_2dfaces(Cells.nodeid, nbelements, faces, cellf)
    elif dim == 3:
        create_3dfaces(Cells.nodeid, nbelements, faces, cellf)
    
    Faces.nodeid, oldTonewIndex = unique(sort(faces), axis=0, return_inverse=True)
    
    nbfaces = len(Faces.nodeid)
    
    Cells.faceid = zeros((nbelements, (dim+1)), dtype=int64)
    create_cell_faceid(nbelements, oldTonewIndex, cellf, Cells.faceid, dim=dim)
    
    ############################################################################
    #creater cells left and right of each face (pycceliser)
    Faces.cellid = -1*ones((nbfaces, 2), dtype=int64)
    create_cellsOfFace(Cells.faceid, nbelements, nbfaces, Faces.cellid, dim=dim)
    ############################################################################
    Cells.cellfid = create_NeighborCellByFace(Cells.faceid, Faces.cellid, nbelements, dim=dim)
    
    ############################################################################
    #create info of faces (pycceliser)
    Faces.name   = zeros(nbfaces, dtype=int64)
    Faces.normal = zeros((nbfaces, 3), dtype=double)
    Faces.mesure = zeros(nbfaces, dtype=double)
    Faces.center = zeros((nbfaces, 3), dtype=double)
    if dim == 2:
        create_info_2dfaces(Faces.cellid, Faces.nodeid, Nodes.name, Nodes.vertex, Cells.center, 
                            nbfaces, Faces.normal, Faces.mesure, Faces.center, Faces.name)
    elif dim == 3:
        create_info_3dfaces(Faces.cellid, Faces.nodeid, Nodes.name, Nodes.vertex, Cells.center, 
                            nbfaces, Faces.normal, Faces.mesure, Faces.center, Faces.name)

    #Number of boundary faces
    Faces.bound = len(Faces.name[Faces.name !=0])
    ############################################################################
    #Create outgoing normal vectors (pycceliser)
    Cells.nf = zeros((nbelements, dim+1, 3), dtype=double)
    create_NormalFacesOfCell(Cells.center, Faces.center, Cells.faceid, Faces.normal, nbelements, Cells.nf, dim=dim)
    ###########################################################################   
    #still numba function
    if dim == 2:
        Faces.oppnodeid = create_2doppNodeOfFaces(Cells.nodeid, Cells.faceid, Faces.nodeid, nbelements, nbfaces)
    elif dim == 3:
        #TODO active if needed
#        Faces.oppnodeid = create_3doppNodeOfFaces(Cells.nodeid, Cells.faceid, Faces.nodeid, nbelements, nbfaces)
        #pycceliser (not yet)
        Faces.nodeid = oriente_3dfacenodeid(Faces.nodeid, Faces.normal, Nodes.vertex)

        
    return 0

def create_2d_halo_structure(file, size):

    for line in file:
        if "endhalosint" in line:
            break
        Halo.halosint.append(int(line))

    for line in file:
        # process the line
        if "halosext" in line:
            continue
        if "centvol" in line:
            break
        Halo.halosext.append([int(x) for x in line.split()])

    for line in file:
        if "centvol" in line:
            continue
        if "globalcelltolocal" in line:
            break
        Halo.centvol.append([float(x) for x in line.split()])

    cmpt = 0
    for line in file:
    #read Global cell to local
        if line == "globalcelltomocal\n":
            continue
        if line == "endglobalcelltolocal\n":
            break
        Cells.globtoloc[int(line)] = cmpt
        Cells.loctoglob.append(int(line))
        cmpt += 1

    cmpt = 0
    for line in file:
        #read Local Node To Global
        if line == "localnodetoglobal\n":
            continue
        if line == "endlocalnodetoglobal\n":
            break
        Nodes.loctoglob.append(int(line))
        cmpt += 1
    for line in file:
        #read LocalToGlobal
        if line == "neigh\n":
            continue
        if line == "endneigh\n":
            break
        Halo.neigh.append([int64(x) for x in line.split()])
    
    for line in file:
        #read LocalToGlobal
        if line == "nodeparts\n":
            continue
        if line == "endnodeparts\n":
            break
        Nodes.nparts.append([int64(x) for x in line.split()])
    
    facenameoldf = Faces.name
    #TODO change adding cell index
    Faces.halofid = zeros(len(Faces.nodeid), dtype=int)
    if size > 1:
        k = 1
        for i in range(len(Halo.halosext)):
            Halo.faces[tuple([Halo.halosext[i][1], Halo.halosext[i][2]])] = k
            Halo.faces[tuple([Halo.halosext[i][2], Halo.halosext[i][3]])] = k + 1
            Halo.faces[tuple([Halo.halosext[i][1], Halo.halosext[i][3]])] = k + 2
            k = k+3

        for i in range(len(Faces.nodeid)):
            n1 = Nodes.loctoglob[Faces.nodeid[i][0]]
            n2 = Nodes.loctoglob[Faces.nodeid[i][1]]
            if  Halo.faces.get(tuple([n1, n2])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10
                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n1, n2])))/3)

            if Halo.faces.get(tuple([n2, n1])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10

                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n2, n1])))/3)
                
        Halo.halosext = asarray(Halo.halosext)
        longueur = 0
        longh = []
        tmp = [[] for i in range(len(Nodes.name))]
        for i in range(len(Nodes.name)):
            if Nodes.name[i] == 10:
                Halo.nodes[i] = Nodes.loctoglob[i]
                arg = where(Halo.halosext[:,1:4] == Nodes.loctoglob[i])
                tmp[i].append(arg[0])
                longueur = max(longueur, len(arg[0]))
                longh.append(len(arg[0]))
            else:
                tmp[i].append(array([-1]))
                longh.append(0)

        Nodes.halonid = [[-1]*longueur for i in range(len(Nodes.name))]
        
        for i in range(len(tmp)):
            for j in range(len(tmp[i][0])):
                Nodes.halonid[i][j] = tmp[i][0][j]
                
        for i in range(len(Nodes.halonid)):
            Nodes.halonid[i].append(longh[i])
                
    if size == 1 : 
        Nodes.halonid = zeros((len(Nodes.name),2), dtype=int)
        Halo.centvol  = zeros((2,2))
        Halo.halosint = zeros((2,2))
        
      # A vérifier !!!!!!
    Faces.ghostcenter = [[] for i in range(len(Faces.name))]
    Nodes.ghostcenter = [[] for i in range(len(Nodes.name))]
    Nodes.haloghostcenter = [[] for i in range(len(Nodes.name))]
    #compute the ghostcenter for each face and each node
    for i in range(len(Faces.name)):
        nod1 = Faces.nodeid[i][1]
        nod2 = Faces.nodeid[i][0]
        if Faces.name[i] == 1 or  Faces.name[i] == 2 or Faces.name[i] ==  3 or Faces.name[i] == 4:
            
            x_1 = Nodes.vertex[nod1]
            x_2 = Nodes.vertex[nod2]
            
            c_left = Faces.cellid[i][0]
            v_1 = Cells.center[c_left]
            gamma = ((v_1[0] - x_2[0])*(x_1[0]-x_2[0]) + (v_1[1]-x_2[1])*(x_1[1]-x_2[1]))/((x_1[0]-x_2[0])**2 + (x_1[1]-x_2[1])**2)
           
            kk = array([gamma * x_1[0] + (1 - gamma) * x_2[0], gamma * x_1[1] + (1 - gamma) * x_2[1]])
            
            v_2 = array([2 * kk[0] + ( -1 * v_1[0]), 2 * kk[1] + ( -1 * v_1[1])])

            Faces.ghostcenter[i] = [v_2[0], v_2[1], gamma]
            Nodes.ghostcenter[nod1].append([v_2[0], v_2[1], Faces.cellid[i][0], facenameoldf[i], i])
            Nodes.ghostcenter[nod2].append([v_2[0], v_2[1], Faces.cellid[i][0], facenameoldf[i], i])
            
        else:
            Faces.ghostcenter[i] = [-1., -1., -1.]
            
    # #update
    # for i in range(len(Nodes.name)):
    #     if len(Nodes.ghostcenter[i]) == 0 :
    #         Nodes.ghostcenter[i].append([0.,0.,-1.])
    #         Nodes.ghostcenter[i].append([0.,0.,-1.])
    #     elif len(Nodes.ghostcenter[i]) == 1 :
    #         Nodes.ghostcenter[i].append([0.,0.,-1.])
           
    Faces.ghostcenter = asarray(Faces.ghostcenter)
    # Nodes.ghostcenter = asarray(Nodes.ghostcenter)
    
    #define halo cells neighbor by nodes
    maxhalonid = 0
    Cells.halonid = [[] for i in range(len(Cells.nodeid))]
    for i in range(len(Cells.nodeid)):
        for j in range(3):
            nod = Cells.nodeid[i][j]
            k = Nodes.halonid[nod][-1]
            Cells.halonid[i].extend(Nodes.halonid[nod][:k])
        Cells.halonid[i] = list(set(Cells.halonid[i]))
        maxhalonid = max(maxhalonid, len(Cells.halonid[i]))
    
    for i in range(len(Cells.nodeid)):
        numb = len(Cells.halonid[i])
        iterator = maxhalonid - len(Cells.halonid[i])
        for k in range(iterator):
             Cells.halonid[i].append(-1)
        Cells.halonid[i].append(numb)
        
    ##########################################################################
    nbelements = len(Cells.volume)
    nbfaces = len(Faces.name)
    nbnodes = len(Nodes.name)
    #Periodic boundary
    #TODO periodic test
    Nodes.periodicid = [[] for i in range(nbnodes)]
    Cells.periodicnid = [[] for i in range(nbelements)]
    Cells.periodicfid = zeros(nbelements)
    Cells.shift = zeros((nbelements, 2))
    
    Nodes.periodicid = asarray(Nodes.periodicid, dtype=int)
    Cells.periodicnid = asarray(Cells.periodicnid, dtype=int)
    Cells.periodicfid = asarray(Cells.periodicfid, dtype=int)
    
    
    #  TODO Periodic boundary (left and right)
    leftb = {}
    rightb = {}
    for i in range(0, nbfaces, 1):
        if Faces.name[i] == 5:
            leftb[float32(Faces.center[i][1])] = Faces.cellid[i][0]
        if Faces.name[i] == 6:
            rightb[float32(Faces.center[i][1])] = Faces.cellid[i][0]
    
    Nodes.periodicid = [[] for i in range(nbnodes)]
    Cells.periodicnid = [[] for i in range(nbelements)]
    Cells.periodicfid = zeros(nbelements)
    Cells.shift = zeros((nbelements, 2))
    
    shiftx = max(Nodes.vertex[:,0])
    shifty = max(Nodes.vertex[:,1])
    
    for i in range(nbfaces):
        if Faces.name[i] == 5:
            
            try:
               
               Faces.cellid[i][1] = rightb[float32(Faces.center[i][1])]
               Cells.periodicfid[Faces.cellid[i][0]] = rightb[float32(Faces.center[i][1])]
               Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[rightb[float32(Faces.center[i][1])]][0:-1])
               for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                   if i != -1:
                       Cells.shift[cell][0] = -1*shiftx

            except KeyError:
              var = 0              
            
            
            
                
                    
            
        elif Faces.name[i] == 6:


            try:
               
               Faces.cellid[i][1] = leftb[float32(Faces.center[i][1])]
               Cells.periodicfid[Faces.cellid[i][0]] = leftb[float32(Faces.center[i][1])]
               Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[leftb[float32(Faces.center[i][1])]][0:-1])
               for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                   if i != -1:
                       Cells.shift[cell][0] = shiftx

            except KeyError:
              var = 0  
    
    leftb = {}
    rightb = {}
    
    for i in range(0, nbnodes, 1):
        for j in range(Nodes.cellid[i][-1]):
            if Nodes.name[i] == 5:
                leftb.setdefault(float32(Nodes.vertex[i][1]), []).append(Nodes.cellid[i][j])
            if Nodes.name[i] == 6:
                rightb.setdefault(float32(Nodes.vertex[i][1]), []).append(Nodes.cellid[i][j])
    
    for i in range(nbnodes):
        if Nodes.name[i] == 5:

            try:
               Nodes.periodicid[i].extend(rightb[float32(Nodes.vertex[i][1])])
               Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]

            except KeyError:
              var = 0 
            
            #Nodes.periodicid[i].extend(rightb[float32(Nodes.vertex[i][1])])
            #Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]
            
        if Nodes.name[i] == 6:

            try:
               Nodes.periodicid[i].extend(leftb[float32(Nodes.vertex[i][1])])
               Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]

            except KeyError:
              var = 0  

            #Nodes.periodicid[i].extend(leftb[float32(Nodes.vertex[i][1])])
            #Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]
    ########################################################################################################
    
    #  TODO Periodic boundary (top and bottom)
    top = {}
    bottom = {}
    for i in range(0, nbfaces, 1):
        if Faces.name[i] == 7:
            top[float32(Faces.center[i][0])] = Faces.cellid[i][0]
        if Faces.name[i] == 8:
            bottom[float32(Faces.center[i][0])] = Faces.cellid[i][0]
    
    for i in range(nbfaces):
        if Faces.name[i] == 7:
        
            try:
               
               Faces.cellid[i][1] = bottom[float32(Faces.center[i][0])]
               Cells.periodicfid[Faces.cellid[i][0]] = bottom[float32(Faces.center[i][0])]
               Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[bottom[float32(Faces.center[i][0])]][0:-1])
               for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                   if i != -1:
                       Cells.shift[cell][1] = shifty

            except KeyError:
              var = 0  
     
        elif Faces.name[i] == 8:
  
            
            try:
               
               Faces.cellid[i][1] = top[float32(Faces.center[i][0])]
               Cells.periodicfid[Faces.cellid[i][0]] = top[float32(Faces.center[i][0])]
               Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[top[float32(Faces.center[i][0])]][0:-1])
               for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                   if i != -1:
                       Cells.shift[cell][1] = -1*shifty

            except KeyError:
              var = 0  

    top = {}
    bottom = {}
    
    # minaxisx = min(Nodes.vertex[:,0]);maxaxisx = max(Nodes.vertex[:,0])
    # minaxisy = min(Nodes.vertex[:,1]);maxaxisy = max(Nodes.vertex[:,1])
    
    for i in range(0, nbnodes, 1):
        for j in range(Nodes.cellid[i][-1]):
            if Nodes.name[i] == 7 :#or (list(Nodes.vertex[i][0:2]) == [minaxisx, maxaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, maxaxisy]):
                top.setdefault(float32(Nodes.vertex[i][0]), []).append(Nodes.cellid[i][j])
            
            if Nodes.name[i] == 8 :#or (list(Nodes.vertex[i][0:2]) == [minaxisx, minaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, minaxisy]):
                bottom.setdefault(float32(Nodes.vertex[i][0]), []).append(Nodes.cellid[i][j])
    
    for i in range(nbnodes):
        if Nodes.name[i] == 7 :#or (list(Nodes.vertex[i][0:2]) == [minaxisx, maxaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, maxaxisy]):
            try:
               Nodes.periodicid[i].extend(bottom[float32(Nodes.vertex[i][0])])

            except KeyError:
              var = 0  
            #Nodes.periodicid[i].extend(bottom[float32(Nodes.vertex[i][0])])
            # Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]
            
        if Nodes.name[i] == 8:# or (list(Nodes.vertex[i][0:2]) == [minaxisx, minaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, minaxisy]):
           # print(Nodes.vertex[i][0])  
            try:
               Nodes.periodicid[i].extend(top[float32(Nodes.vertex[i][0])])

            except KeyError:
              var = 0  
            #Nodes.periodicid[i].extend(top[float32(Nodes.vertex[i][0])])
            
            # Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]
    
    
    maxperiodiccell = 0
    for i in range(nbelements):
        Cells.periodicnid[i] = unique(Cells.periodicnid[i])
        maxperiodiccell = max(maxperiodiccell, len(Cells.periodicnid[i]))
    
    for i in range(nbelements):
        iterator  = maxperiodiccell - len(Cells.periodicnid[i])
        for j in range(iterator):
            Cells.periodicnid[i] = append(Cells.periodicnid[i], -1)
        if Cells.periodicfid[i] == 0:
            Cells.periodicfid[i] = -1
     
        
    maxperiodicnode = 0
    for i in range(nbnodes):
        # Nodes.periodicid[i] = np.unique(Nodes.periodicid[i])
        maxperiodicnode = max(maxperiodicnode, len(Nodes.periodicid[i]))
    
    for i in range(nbnodes):
        iterator  = maxperiodicnode - len(Nodes.periodicid[i])
        for j in range(iterator):
            Nodes.periodicid[i] = append(Nodes.periodicid[i], -1) 
        
            
    Nodes.periodicid = asarray(Nodes.periodicid, dtype=int)
    Cells.periodicnid = asarray(Cells.periodicnid, dtype=int)
    Cells.periodicfid = asarray(Cells.periodicfid, dtype=int)
        
      
    if size == 1 :
        Cells.halonid = zeros((len(Cells.nodeid),2), dtype=int)
        Halo.halosext = zeros((2,2), dtype=int)
        for i in range(len(Nodes.vertex)):
            Nodes.loctoglob.append(i)
        
        for i in range(len(Cells.nodeid)):
             Cells.loctoglob.append(i)
        
    cells = Cells(Cells.nodeid,  Cells.faceid, Cells.center,
                  Cells.volume,  Cells.cellfid, Cells.cellnid, 
                  Cells.nf, Cells.globtoloc, asarray(Cells.loctoglob), asarray(Cells.halonid), 
                  asarray(Cells.tc), Cells.shift, Cells.periodicfid, Cells.periodicnid)
    
    nodes = Nodes(Nodes.vertex, Nodes.name, Nodes.cellid, 
                  Nodes.ghostcenter, Nodes.loctoglob, asarray(Nodes.halonid),
                  Nodes.nparts, Nodes.haloghostcenter, Nodes.periodicid)

    faces = Faces(Faces.nodeid, Faces.cellid, Faces.name,
                  Faces.normal, Faces.mesure, Faces.center, Faces.bound, 
                  Faces.ghostcenter, Faces.oppnodeid, asarray(Faces.halofid)) 
    
    halos = Halo(asarray(Halo.halosint), asarray(Halo.halosext), asarray(Halo.centvol),
                 asarray(Halo.neigh), Halo.faces, Halo.nodes)
    
    return cells, faces, nodes, halos


def create_3d_halo_structure(file, size):
    
    for line in file:
        if "endhalosint" in line:
            break
        Halo.halosint.append(int(line))#.append(int(line))#for x in line.split()])

    for line in file:
        # process the line
        if "halosext" in line:
            continue
        if "centvol" in line:
            break
        Halo.halosext.append([int(x) for x in line.split()])

    for line in file:
        if "centvol" in line:
            continue
        if "globalcelltolocal" in line:
            break
        Halo.centvol.append([float(x) for x in line.split()])

    cmpt = 0
    for line in file:
    #read Global cell to local
        if line == "globalcelltomocal\n":
            continue
        if line == "endglobalcelltolocal\n":
            break
        Cells.globtoloc[int(line)] = cmpt
        Cells.loctoglob.append(int(line))
        cmpt += 1

    cmpt = 0
    for line in file:
        #read Local Node To Global
        if line == "localnodetoglobal\n":
            continue
        if line == "endlocalnodetoglobal\n":
            break
        Nodes.loctoglob.append(int(line))
        # cmpt += 1
    for line in file:
        #read LocalToGlobal
        if line == "neigh\n":
            continue
        if line == "endneigh\n":
            break
        Halo.neigh.append([int(x) for x in line.split()])
    for line in file:
        #read LocalToGlobal
        if line == "nodeparts\n":
            continue
        if line == "endnodeparts\n":
            break
        Nodes.nparts.append([int64(x) for x in line.split()])
        
    Faces.halofid = zeros(len(Faces.nodeid), dtype=int)
    
    facenameoldf = Faces.name
    #TODO change adding cell index
    if size > 1:
        k = 1
        for i in range(len(Halo.halosext)):
            Halo.faces[tuple([Halo.halosext[i][1], Halo.halosext[i][2],  Halo.halosext[i][3] ])] = k
            Halo.faces[tuple([Halo.halosext[i][3], Halo.halosext[i][4],  Halo.halosext[i][1] ])] = k + 1
            Halo.faces[tuple([Halo.halosext[i][1], Halo.halosext[i][2],  Halo.halosext[i][4] ])] = k + 2
            Halo.faces[tuple([Halo.halosext[i][4], Halo.halosext[i][2],  Halo.halosext[i][3] ])] = k + 3
            
            k = k+4

        for i in range(len(Faces.nodeid)):
            n1 = Nodes.loctoglob[Faces.nodeid[i][0]]
            n2 = Nodes.loctoglob[Faces.nodeid[i][1]]
            n3 = Nodes.loctoglob[Faces.nodeid[i][2]]
            
            if  Halo.faces.get(tuple([n1, n2, n3])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10
                Nodes.name[Faces.nodeid[i][2]] = 10
               
                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n1, n2, n3])))/4)
            
            if Halo.faces.get(tuple([n1, n3, n2])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10
                Nodes.name[Faces.nodeid[i][2]] = 10
                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n1, n3, n2])))/4)

            if Halo.faces.get(tuple([n2, n1, n3])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10
                Nodes.name[Faces.nodeid[i][2]] = 10
                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n2, n1, n3])))/4)
                
            if Halo.faces.get(tuple([n2, n3, n1])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10
                Nodes.name[Faces.nodeid[i][2]] = 10
                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n2, n3, n1])))/4)
                
            if Halo.faces.get(tuple([n3, n1, n2])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10
                Nodes.name[Faces.nodeid[i][2]] = 10
                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n3, n1, n2])))/4)
            
            if Halo.faces.get(tuple([n3, n2, n1])):

                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.name[i] = 10
                Nodes.name[Faces.nodeid[i][0]] = 10
                Nodes.name[Faces.nodeid[i][1]] = 10
                Nodes.name[Faces.nodeid[i][2]] = 10
                Faces.halofid[i] = int((-1+Halo.faces.get(tuple([n3, n2, n1])))/4)
        
        #TODO verify halosext
        Halo.halosext = asarray(Halo.halosext)
        longueur = 0
        longh = []
        tmp = [[] for i in range(len(Nodes.name))]
        for i in range(len(Nodes.name)):
            if Nodes.name[i] == 10:
                Halo.nodes[i] = Nodes.loctoglob[i]
                arg = where(Halo.halosext[:,1:5] == Nodes.loctoglob[i])
                tmp[i].append(arg[0])
                longueur = max(longueur, len(arg[0]))
                longh.append(len(arg[0]))
            else:
                tmp[i].append(array([-1]))
                longh.append(0)

        Nodes.halonid = [[-1]*longueur for i in range(len(Nodes.name))]
        
        for i in range(len(tmp)):
            for j in range(len(tmp[i][0])):
                Nodes.halonid[i][j] = tmp[i][0][j]
                
        for i in range(len(Nodes.halonid)):
            Nodes.halonid[i].append(longh[i])
                
    if size == 1 : 
        Nodes.halonid = zeros((len(Nodes.name),2), dtype=int)
        Halo.centvol  = zeros((2,2))
        Halo.halosint = zeros((2,2))
        for i in range(len(Cells.nodeid)):
             Cells.loctoglob.append(i)
        
        
    Faces.ghostcenter = zeros((len(Faces.name), 4))#[[] for i in range(len(Faces.name))]
    Nodes.ghostcenter = [[] for i in range(len(Nodes.name))]
    Nodes.haloghostcenter = [[] for i in range(len(Nodes.name))]
     
    kk = zeros(3)
    #TODO ghost center à verifier
    #compute the ghostcenter for each face and each node
    for i in range(len(Faces.name)):
        
        if Faces.name[i] != 0 and Faces.name[i] != 10:
            nod1 = Faces.nodeid[i][1]
            nod2 = Faces.nodeid[i][0]
            nod3 = Faces.nodeid[i][2]
            
            n = Faces.normal[i]/Faces.mesure[i]
            
            c_left = Faces.cellid[i][0]
            v_1 = Cells.center[c_left]
            u = Faces.center[i][:] - v_1[:]
            gamma = dot(u, n)
            
            kk[0] = v_1[0] + gamma*n[0]
            kk[1] = v_1[1] + gamma*n[1]
            kk[2] = v_1[2] + gamma*n[2]
            
            v_2 = array([2 * kk[0] + ( -1 * v_1[0]), 2 * kk[1] + ( -1 * v_1[1]), 2 * kk[2] + ( -1 * v_1[2])])
            
            Faces.ghostcenter[i] = [v_2[0], v_2[1], v_2[2], gamma]
            Nodes.ghostcenter[nod1].append([v_2[0], v_2[1], v_2[2], Faces.cellid[i][0], facenameoldf[i], i])
            Nodes.ghostcenter[nod2].append([v_2[0], v_2[1], v_2[2], Faces.cellid[i][0], facenameoldf[i], i])
            Nodes.ghostcenter[nod3].append([v_2[0], v_2[1], v_2[2], Faces.cellid[i][0], facenameoldf[i], i])
            
        else:
            Faces.ghostcenter[i] = [0.,0.,0., -1.]
        
    # maxGhostCell = 0
    # for i in range(len(Nodes.name)):
    #     maxGhostCell = max(maxGhostCell, len(Nodes.ghostcenter[i]))
    
    # for i in range(len(Nodes.name)):
    #     iterator = maxGhostCell - len(Nodes.ghostcenter[i])
    #     for k in range(iterator):
    #         Nodes.ghostcenter[i].append([-1., -1., -1., -1.])
    
    #     if len(Nodes.ghostcenter[i]) == 0 :
    #             Nodes.ghostcenter[i].append([-1, -1., -1.,-1.])
    
    # Nodes.ghostcenter = asarray(Nodes.ghostcenter)
    #define halo cells neighbor by nodes
    maxhalonid = 0
    Cells.halonid = [[] for i in range(len(Cells.nodeid))]
    for i in range(len(Cells.nodeid)):
        for j in range(4):
            nod = Cells.nodeid[i][j]
            k = Nodes.halonid[nod][-1]
            Cells.halonid[i].extend(Nodes.halonid[nod][:k])
        Cells.halonid[i] = list(set(Cells.halonid[i]))
        maxhalonid = max(maxhalonid, len(Cells.halonid[i]))
    
    for i in range(len(Cells.nodeid)):
        numb = len(Cells.halonid[i])
        iterator = maxhalonid - len(Cells.halonid[i])
        for k in range(iterator):
             Cells.halonid[i].append(-1)
        Cells.halonid[i].append(numb)
        
        
    ##########################################################################
    nbelements = len(Cells.volume)
    nbfaces = len(Faces.name)
    nbnodes = len(Nodes.name)
    #Periodic boundary
    #TODO periodic test
    Nodes.periodicid = [[] for i in range(nbnodes)]
    Cells.periodicnid = [[] for i in range(nbelements)]
    Cells.periodicfid = zeros(nbelements)
    Cells.shift = zeros((nbelements, 2))
    
    Nodes.periodicid = asarray(Nodes.periodicid, dtype=int)
    Cells.periodicnid = asarray(Cells.periodicnid, dtype=int)
    Cells.periodicfid = asarray(Cells.periodicfid, dtype=int)
    
    
    #  TODO Periodic boundary (left and right)
    leftb = {}
    rightb = {}
    for i in range(0, nbfaces, 1):
        if Faces.name[i] == 5:
            leftb[float32(Faces.center[i][1])] = Faces.cellid[i][0]
        if Faces.name[i] == 6:
            rightb[float32(Faces.center[i][1])] = Faces.cellid[i][0]
    
    Nodes.periodicid = [[] for i in range(nbnodes)]
    Cells.periodicnid = [[] for i in range(nbelements)]
    Cells.periodicfid = zeros(nbelements)
    Cells.shift = zeros((nbelements, 2))
    
    shiftx = max(Nodes.vertex[:,0])
    shifty = max(Nodes.vertex[:,1])
    
    for i in range(nbfaces):
        if Faces.name[i] == 5:
            
            Faces.cellid[i][1] = rightb[float32(Faces.center[i][1])]
            Cells.periodicfid[Faces.cellid[i][0]] = rightb[float32(Faces.center[i][1])]
            Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[rightb[float32(Faces.center[i][1])]][0:-1])
            for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                if i != -1:
                    Cells.shift[cell][0] = -1*shiftx
            
        elif Faces.name[i] == 6:
            Faces.cellid[i][1] = leftb[float32(Faces.center[i][1])]
            Cells.periodicfid[Faces.cellid[i][0]] = leftb[float32(Faces.center[i][1])]
            Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[leftb[float32(Faces.center[i][1])]][0:-1])
            for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                if i != -1:
                    Cells.shift[cell][0] = shiftx
    
    leftb = {}
    rightb = {}
    
    for i in range(0, nbnodes, 1):
        for j in range(Nodes.cellid[i][-1]):
            if Nodes.name[i] == 5:
                leftb.setdefault(float32(Nodes.vertex[i][1]), []).append(Nodes.cellid[i][j])
            if Nodes.name[i] == 6:
                rightb.setdefault(float32(Nodes.vertex[i][1]), []).append(Nodes.cellid[i][j])
    
    for i in range(nbnodes):
        if Nodes.name[i] == 5:
            Nodes.periodicid[i].extend(rightb[float32(Nodes.vertex[i][1])])
           
            
        if Nodes.name[i] == 6:
            Nodes.periodicid[i].extend(leftb[float32(Nodes.vertex[i][1])])
    ########################################################################################################
    
    #  TODO Periodic boundary (top and bottom)
    top = {}
    bottom = {}
    for i in range(0, nbfaces, 1):
        if Faces.name[i] == 7:
            top[float32(Faces.center[i][0])] = Faces.cellid[i][0]
        if Faces.name[i] == 8:
            bottom[float32(Faces.center[i][0])] = Faces.cellid[i][0]
    
    for i in range(nbfaces):
        if Faces.name[i] == 7:
            Faces.cellid[i][1] = bottom[float32(Faces.center[i][0])]
            Cells.periodicfid[Faces.cellid[i][0]] = bottom[float32(Faces.center[i][0])]
            Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[bottom[float32(Faces.center[i][0])]][0:-1])
            for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                if i != -1:
                    Cells.shift[cell][1] = shifty
            
        elif Faces.name[i] == 8:
            Faces.cellid[i][1] = top[float32(Faces.center[i][0])]
            Cells.periodicfid[Faces.cellid[i][0]] = top[float32(Faces.center[i][0])]
            Cells.periodicnid[Faces.cellid[i][0]].extend(Cells.cellnid[top[float32(Faces.center[i][0])]][0:-1])
            for cell in Cells.periodicnid[Faces.cellid[i][0]]:
                if i != -1:
                    Cells.shift[cell][1] = -1*shifty
    
    top = {}
    bottom = {}
    
    # minaxisx = min(Nodes.vertex[:,0]);maxaxisx = max(Nodes.vertex[:,0])
    # minaxisy = min(Nodes.vertex[:,1]);maxaxisy = max(Nodes.vertex[:,1])
    
    for i in range(0, nbnodes, 1):
        for j in range(Nodes.cellid[i][-1]):
            if Nodes.name[i] == 7 :#or (list(Nodes.vertex[i][0:2]) == [minaxisx, maxaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, maxaxisy]):
                top.setdefault(float32(Nodes.vertex[i][0]), []).append(Nodes.cellid[i][j])
            
            if Nodes.name[i] == 8 :#or (list(Nodes.vertex[i][0:2]) == [minaxisx, minaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, minaxisy]):
                bottom.setdefault(float32(Nodes.vertex[i][0]), []).append(Nodes.cellid[i][j])
    
    for i in range(nbnodes):
        if Nodes.name[i] == 7 :#or (list(Nodes.vertex[i][0:2]) == [minaxisx, maxaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, maxaxisy]):
            Nodes.periodicid[i].extend(bottom[float32(Nodes.vertex[i][0])])
            # Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]
            
        if Nodes.name[i] == 8:# or (list(Nodes.vertex[i][0:2]) == [minaxisx, minaxisy]) or \
                # (list(Nodes.vertex[i][0:2]) == [maxaxisx, minaxisy]):
            Nodes.periodicid[i].extend(top[float32(Nodes.vertex[i][0])])
            # Nodes.ghostcenter[i] = [[-1., -1., -1., -1., -1.]]
    
    
    maxperiodiccell = 0
    for i in range(nbelements):
        Cells.periodicnid[i] = unique(Cells.periodicnid[i])
        maxperiodiccell = max(maxperiodiccell, len(Cells.periodicnid[i]))
    
    for i in range(nbelements):
        iterator  = maxperiodiccell - len(Cells.periodicnid[i])
        for j in range(iterator):
            Cells.periodicnid[i] = append(Cells.periodicnid[i], -1)
        if Cells.periodicfid[i] == 0:
            Cells.periodicfid[i] = -1
     
        
    maxperiodicnode = 0
    for i in range(nbnodes):
        # Nodes.periodicid[i] = np.unique(Nodes.periodicid[i])
        maxperiodicnode = max(maxperiodicnode, len(Nodes.periodicid[i]))
    
    for i in range(nbnodes):
        iterator  = maxperiodicnode - len(Nodes.periodicid[i])
        for j in range(iterator):
            Nodes.periodicid[i] = append(Nodes.periodicid[i], -1) 
        
            
    Nodes.periodicid = asarray(Nodes.periodicid, dtype=int)
    Cells.periodicnid = asarray(Cells.periodicnid, dtype=int)
    Cells.periodicfid = asarray(Cells.periodicfid, dtype=int)

    
        
    if size == 1 :
        Cells.halonid = zeros((len(Cells.nodeid),2), dtype=int)
        Halo.halosext = zeros((2,2), dtype=int)
        for i in range(len(Nodes.vertex)):
            Nodes.loctoglob.append(i)
        
    cells = Cells(Cells.nodeid,  Cells.faceid, Cells.center,
                  Cells.volume,  Cells.cellfid, Cells.cellnid, 
                  Cells.nf, Cells.globtoloc, asarray(Cells.loctoglob), 
                  asarray(Cells.halonid), asarray(Cells.tc), Cells.shift, Cells.periodicfid, Cells.periodicnid)
    
    nodes = Nodes(Nodes.vertex, Nodes.name, Nodes.cellid, 
                  Nodes.ghostcenter, asarray(Nodes.loctoglob), asarray(Nodes.halonid),
                  Nodes.nparts, Nodes.haloghostcenter, Nodes.periodicid)

    faces = Faces(Faces.nodeid, Faces.cellid, Faces.name,
                  Faces.normal, Faces.mesure, Faces.center, Faces.bound, 
                  Faces.ghostcenter, Faces.oppnodeid, asarray(Faces.halofid)) 
    
    halos = Halo(asarray(Halo.halosint), asarray(Halo.halosext), asarray(Halo.centvol),
                 asarray(Halo.neigh), Halo.faces, Halo.nodes)
    
    return cells, faces, nodes, halos

from mpi4py import MPI

#COMM = MPI.COMM_WORLD
#SIZE = COMM.Get_size()
#RANK = COMM.Get_rank()

def generate_structure(dim, size):
    
    import os

    MESH_DIR = "meshes"+str(size)+"PROC"
    filename = os.path.join(MESH_DIR, 'mesh'+str(RANK)+'.txt')
    txt_file = open(filename)
    
    start = timeit.default_timer()
    CreateStructure(txt_file, dim=dim, size=size)
    
    if dim == 2:
        cells, faces, nodes, halos = create_2d_halo_structure(txt_file, size)
    elif dim == 3:
        cells, faces, nodes, halos = create_3d_halo_structure(txt_file, size)
    stop = timeit.default_timer()
    if RANK==0:
        print("CPU time for creating "+str(dim)+"d structure ", stop-start)
   
    grid = {}

    grid["cells"] = cells
    grid["nodes"] = nodes
    grid["faces"] = faces
    grid["halos"] = halos

    txt_file.close()

    return grid
