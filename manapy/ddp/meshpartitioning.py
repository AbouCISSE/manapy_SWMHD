#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:08:27 2021

@author: kissami
"""


# coding: utf-8
import os
import timeit
import meshio
import numpy as np

__all__ = ['readmesh']
from collections import OrderedDict
from mgmetis import metis
from mgmetis.enums import OPTION
from numba import njit

@njit
def create_npart_cpart(cell_nodeid:'int[:,:]', npart:'int[:]', epart:'int[:]', nbnodes:'int', nbelements:'int',
                       size:'int', dim:'int'):
    npart       = [ [i ] for i in npart ]
    cpart       = [ [i ] for i in epart ]
    neighsub    = [[i for i in range(0)]  for i in range(size)]
    cell_part   = [[i for i in range(0)]  for i in range(size)]
    node_part   = [[i for i in range(0)]  for i in range(size)]
    halo_cellid = [[i for i in range(0)]  for i in range(size)]

    for i in range(nbelements):
        for j in range(dim+1):
            k = cell_nodeid[i][j]
            if epart[i] not in npart[k]:
                npart[k].append(epart[i])
#            if k not in node_part[epart[i]]:
            node_part[epart[i]].append(k)
        cell_part[epart[i]].append(i)
    
    
    for i in range(nbelements):
        for j in range(dim+1):
            for k in range(len(npart[cell_nodeid[i][j]])):
                if npart[cell_nodeid[i][j]][k] not in cpart[i]:
                    cpart[i].append(npart[cell_nodeid[i][j]][k])
                    
    for i in range(nbnodes):
        if len(npart[i]) > 1:
            for j in range(len(npart[i])):
                neighsub[npart[i][j]].extend(npart[i])
    
    for i in range(size):
        for j in cell_part[i]:
            a = 0.
            for k in range(dim+1):
                a += len(npart[cell_nodeid[j][k]])
                if a > dim+1:
                    halo_cellid[i].append(j)
                
    return npart, cpart, neighsub, node_part, cell_part, halo_cellid

def unique_func(array):
        uniq, index = np.unique(array, return_index=True)
        return uniq[index.argsort()]
def wedge_3d(u, v):
    wedge = np.zeros(3)
    
    wedge[0] = u[1]*v[2] - u[2]*v[1]
    wedge[1] = u[2]*v[0] - u[0]*v[2]
    wedge[2] = u[0]*v[1] - u[1]*v[0]
    
    return wedge

def dot_vec3(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def det_vec_3d(u, v, w):
    return dot_vec3(u, wedge_3d(v, w))

def readmesh(filename, dim, size, periodic):

    if dim == 2 :
        typeOfCells = "triangle"
        typeOfFaces = "line"
    else:
        typeOfCells = "tetra"
        typeOfFaces = "triangle"

    def load_gmsh_mesh(filename):
        #mesh = meshio.gmsh.read(filename)
        mesh = meshio.read(filename)
        return mesh

    def create_cell_nodeid(mesh):
        cell_nodeid = []
        
        if type(mesh.cells) == dict:
            cell_nodeid = mesh.cells[typeOfCells]

        elif type(mesh.cells) == list:
            cell_nodeid = mesh.cells[1].data

        for i in range(len(cell_nodeid)):
            cell_nodeid[i].sort()
            
        return cell_nodeid

    def define_ghost_node(mesh, nodes):
        ghost_nodes = [0]*len(nodes)

        if type(mesh.cells) == dict:
            for i, j in mesh.cell_data.items():
                if i == typeOfFaces:
                    ghost = j.get('gmsh:physical')

            for i, j in mesh.cells.items():
                if i == typeOfFaces:
                    for k in range(len(j)):
                        for index in range(dim):
                            if ghost[k] == 1 or ghost[k] == 2 or ghost[k] == 5 or ghost[k] == 6 :
                                ghost_nodes[j[k][index]] = int(ghost[k])
            
            
            for i, j in mesh.cells.items():
                if i == typeOfFaces:
                    for k in range(len(j)):
                        for index in range(dim):
                            if ghost_nodes[j[k][index]] == 0:
                                ghost_nodes[j[k][index]] = int(ghost[k])

        elif type(mesh.cells) == list:
            ghost = mesh.cell_data['gmsh:physical'][0]
            for i in range(len(mesh.cells[0].data)):
                for j in range(dim):
                    if ghost[k] == 1 or ghost[k] == 2 or ghost[k] == 5 or ghost[k] == 6 :
                        ghost_nodes[mesh.cells[0].data[i][j]] = int(ghost[i])

            for i in range(len(mesh.cells[0].data)):
                for j in range(dim):
                    if ghost_nodes[j[k][index]] == 0:
                        ghost_nodes[mesh.cells[0].data[i][j]] = int(ghost[i])
                        
        if periodic[0] == 1:
            for i in range(len(ghost_nodes)):
                if ghost_nodes[i] == 1:
                    ghost_nodes[i] = 5
                elif ghost_nodes[i] == 2:
                    ghost_nodes[i] = 6
                    
        if periodic[1] == 1:
            for i in range(len(ghost_nodes)):
                if ghost_nodes[i] == 3:
                    ghost_nodes[i] = 7
                elif ghost_nodes[i] == 4:
                    ghost_nodes[i] = 8

        return ghost_nodes

    def create_nodes(mesh):
        nodes = []
        nodes = mesh.points
        return nodes
    print("Starting ....")
    start = timeit.default_timer()

    #load mesh
    mesh = load_gmsh_mesh(filename)
    #coordinates x, y of each node
    nodes = create_nodes(mesh)
    #nodes of each cell
    cell_nodeid = create_cell_nodeid(mesh)

    ghost_nodes = define_ghost_node(mesh, nodes)
    
    nbelements = len(cell_nodeid)
    nbnodes = len(nodes)

    print("Number of Cells : ", nbelements)
    print("Number of Nodes : ", nbnodes)
    
    MESH_DIR = "meshes"+str(size)+"PROC"
    if not os.path.exists(MESH_DIR):
        os.mkdir(MESH_DIR)
    
    if size == 1:
        for i in range(size):
            # THe mesh file
            filename = os.path.join(MESH_DIR, "mesh"+str(i)+".txt")
            if os.path.exists(filename):
                os.remove(filename)
        
            if os.path.exists(filename):
                os.remove(filename)
    
        with open(filename, "a") as text_file:
            text_file.write("elements\n")
            np.savetxt(text_file, cell_nodeid, fmt='%u')
            text_file.write("endelements\n")
            text_file.write("nodes\n")
            for i in range(len(nodes)):
                for j in range(3):
                    text_file.write(str(nodes[i][j])+str(" "))
                text_file.write(str(ghost_nodes[i]))
                text_file.write("\n")
            text_file.write("endnodes\n")
    
        stop = timeit.default_timer()
        print('CPU time for generating '+str(dim)+'d mesh : ', stop - start)
    else :
        startmetis = timeit.default_timer()
        #Partitioning mesh
        opts = metis.get_default_options()
        opts[OPTION.MINCONN] = 1
        opts[OPTION.CONTIG] = 1
        #options = opts
	
        objval, epart, npart = metis.part_mesh_dual(size, cell_nodeid, nv=len(mesh.points))
	        
    
        stopmetis = timeit.default_timer()
        print("METIS partitionning in ", size, "partitions", stopmetis - startmetis)  
        globnodetoloc = OrderedDict()
        locnodetoglob = OrderedDict()
        globcelltoloc = OrderedDict()
        
        #TODO improve function 
        npart, cpart, neighsub, node_part, cell_part, halo_cellid = create_npart_cpart(cell_nodeid, npart, epart, 
                                                                                       nbnodes, nbelements, size, dim)
        
        
        for i in range(size): 
            neighsub[i]  = np.unique(neighsub[i])
            cell_part[i] = np.unique(cell_part[i])
            neighsub[i]  = neighsub[i][neighsub[i]!=i]
            node_part[i] = unique_func(node_part[i])
            halo_cellid[i]  = np.unique(halo_cellid[i])
            
        for i in range(size):
            for j in range(len(cell_part[i])):
                globcelltoloc[i, j] = cell_part[i][j]
            for j in range(len(node_part[i])):
                globnodetoloc[i, node_part[i][j]] = j
                locnodetoglob[i, j] = node_part[i][j]
        
        cmpt = 0
        tc = np.zeros(nbelements, dtype=np.int32)
        for i in range(size):
            for j in range(len(cell_part[i])):
                tc[cmpt] = cell_part[i][j]
                cmpt += 1
        stopstruc = timeit.default_timer()
        print("Create local structure for each proc", stopstruc - stopmetis)
        
        haloint = OrderedDict()
        haloext = OrderedDict()
        for i in range(size):
            for cell in halo_cellid[i]:
                for k in range(len(cpart[cell])):
                    if i != cpart[cell][k]:
                        haloint.setdefault((i, cpart[cell][k]), []).append(cell)
                        haloext.setdefault((cpart[cell][k], i), []).append(cell)
        
        
        centvol = [[] for i in range(size)]
        if dim == 2:
            for i in range(size):
                for j in range(len(neighsub[i])):
                    for k in range(len(haloext[(i, neighsub[i][j])])):
                        s_1 = cell_nodeid[haloext[(i, neighsub[i][j])][k]][0]
                        s_2 = cell_nodeid[haloext[(i, neighsub[i][j])][k]][1]
                        s_3 = cell_nodeid[haloext[(i, neighsub[i][j])][k]][2]
            
                        x_1 = nodes[s_1][0]
                        y_1 = nodes[s_1][1]
                        x_2 = nodes[s_2][0]
                        y_2 = nodes[s_2][1]
                        x_3 = nodes[s_3][0]
                        y_3 = nodes[s_3][1]
            
                        centvol[i].append([1./3 * (x_1 + x_2 + x_3), 1./3*(y_1 + y_2 + y_3), 0.,
                                           (1./2) * abs((x_1-x_2)*(y_1-y_3)-(x_1-x_3)*(y_1-y_2))])
    
        if dim == 3:
            for i in range(size):
                for j in range(len(neighsub[i])):
                    for k in range(len(haloext[(i, neighsub[i][j])])):
                        
                        s_1 = cell_nodeid[haloext[(i, neighsub[i][j])][k]][0]
                        s_2 = cell_nodeid[haloext[(i, neighsub[i][j])][k]][1]
                        s_3 = cell_nodeid[haloext[(i, neighsub[i][j])][k]][2]
                        s_4 = cell_nodeid[haloext[(i, neighsub[i][j])][k]][3]
                        
                        a = np.asarray(nodes[s_1])
                        b = np.asarray(nodes[s_2])
                        c = np.asarray(nodes[s_3])
                        d = np.asarray(nodes[s_4])
            
                        x_1 = nodes[s_1][0]; y_1 = nodes[s_1][1]; z_1 = nodes[s_1][2]
                        x_2 = nodes[s_2][0]; y_2 = nodes[s_2][1]; z_2 = nodes[s_2][2] 
                        x_3 = nodes[s_3][0]; y_3 = nodes[s_3][1]; z_3 = nodes[s_3][2] 
                        x_4 = nodes[s_4][0]; y_4 = nodes[s_4][1]; z_4 = nodes[s_4][2]
            
                        centvol[i].append([1./4 * (x_1 + x_2 + x_3 + x_4), 1./4*(y_1 + y_2 + y_3 + y_4), 
                                           1./4*(z_1 + z_2 + z_3 + z_4), 
                                           1./6*np.fabs(det_vec_3d(b-a, c-a, d-a))])
        
        stophalo = timeit.default_timer()
        print("Creating halo structure", stophalo - stopstruc)
                
        for i in range(size):
            # THe mesh file
            filename = os.path.join(MESH_DIR, "mesh"+str(i)+".txt")
            if os.path.exists(filename):
                os.remove(filename)
        filename = os.path.join(MESH_DIR, "mesh0.txt")
        with open(filename, "a") as text_file:
            text_file.write("GtoL\n")
            for i in range(len(tc)):
                text_file.write(str(tc[i]))
                text_file.write("\n")
            text_file.write("endGtoL\n")
        for i in range(size):
            filename = os.path.join(MESH_DIR, "mesh"+str(i)+".txt")
            with open(filename, "a") as text_file:
                text_file.write("elements\n")
                for j in range(len(cell_part[i])):
                    for k in range(dim+1):
                        text_file.write(str(globnodetoloc[i, cell_nodeid[cell_part[i][j]][k]]))
                        text_file.write(" ")
                    text_file.write("\n")
                text_file.write("endelements\n")
                text_file.write("nodes\n")
                for j in range(len(node_part[i])):
                    text_file.write(str(nodes[node_part[i][j]][0])+" "+str(nodes[node_part[i][j]][1])+ " "+
                                    str(nodes[node_part[i][j]][2])+" "+str(ghost_nodes[node_part[i][j]]))
                    text_file.write("\n")
                text_file.write("endnodes\n")
                text_file.write("halosint\n")
                for j in range(len(neighsub[i])):
                    for k in range(len(haloint[(i, neighsub[i][j])])):
                        text_file.write(str(haloint[(i, neighsub[i][j])][k]))
                        text_file.write("\n")
                text_file.write("endhalosint\n")
                text_file.write("halosext\n")
                for j in range(len(neighsub[i])):
                    for k in range(len(haloext[(i, neighsub[i][j])])):
                        text_file.write(str(haloext[(i, neighsub[i][j])][k])+" " )
                        for m in range(dim+1):
                            text_file.write(str(cell_nodeid[haloext[(i, neighsub[i][j])][k]][m]))
                            text_file.write(" ")
                        text_file.write("\n")
                text_file.write("endhalosext\n")
                text_file.write("centvol\n")
                np.savetxt(text_file, centvol[i])
                text_file.write("endcentvol\n")
                text_file.write("globalcelltolocal\n")
                for j in range(len(cell_part[i])):
                    text_file.write(str(globcelltoloc[i, j]))
                    text_file.write("\n")
                text_file.write("endglobalcelltolocal\n")
                text_file.write("localnodetoglobal\n")
                for j in range(len(node_part[i])):
                    text_file.write(str(locnodetoglob[i, j]))
                    text_file.write("\n")
                text_file.write("endlocalnodetoglobal\n")
                text_file.write("neigh\n")
                for j in range(len(neighsub[i])):
                    text_file.write(str(neighsub[i][j])+ " ")
                text_file.write("\n")
                for j in neighsub[i]:
                    text_file.write(str(len(haloint[(i, j)]))+ " ")
                text_file.write("\n")
                text_file.write("endneigh\n")   
                text_file.write("nodeparts\n")
                for j in range(len(npart)):
                    for k in range(len(npart[j])):
                        text_file.write(str(npart[j][k])+ " ")
                    text_file.write("\n")
                text_file.write("endnodeparts\n")   
                
        stopfile = timeit.default_timer()
        print("save structures in files", stopfile - stophalo)
        
        stop = timeit.default_timer()
        print('CPU time for generating '+str(dim)+'d meshes : ', stop - start)
