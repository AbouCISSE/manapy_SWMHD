#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 09:41:06 2021

@author: kissami
"""
import numpy as np
from mpi4py import MPI
from numba import njit
from ctypes import c_double, POINTER, c_void_p, CDLL, c_int
import os

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

#os.system("mpicc -fPIC -shared -o comm.so comm.c")
try:
    COMMFILE = os.environ['COMMFILE']
except:
    BASE_DIR = os.path.dirname(os.path.realpath(__file__))
    COMMFILE = os.path.join(BASE_DIR, '.')
filename = os.path.join(COMMFILE, "comm.so")

commFunc = CDLL(filename)
commFunc.comm_neigh.argtypes = [POINTER(c_double),POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_int), 
                                POINTER(c_int) , c_void_p]

@njit
def define_halosend(w_c:'float[:]', w_halosend:'float[:]', indsend:'int[:]', nbproc:'int'):
    if nbproc > 1:
        w_halosend[:] = w_c[indsend[:]]

def create_mpi_graph(neighbors):
    topo = COMM.Create_dist_graph_adjacent(neighbors, neighbors,
                                           sourceweights=None, destweights=None)
    return topo
    

def all_to_all(w_halosend, taille, scount, sdepl, rcount, rdepl, w_halorecv, comm_ptr):
    
    address = MPI._addressof(comm_ptr)
    comm_ptr = c_void_p(address)
    commFunc.comm_neigh(w_halosend.ctypes.data_as(POINTER(c_double)), w_halorecv.ctypes.data_as(POINTER(c_double)), 
                           sdepl, scount,  rdepl, rcount, comm_ptr)


def prepare_comm(cells, halos):

    scount = np.zeros(SIZE, dtype=np.intc)
    sdepl = np.zeros(SIZE, dtype=np.intc)
    rcount = np.zeros(SIZE, dtype=np.intc)
    rdepl = np.zeros(SIZE, dtype=np.intc)
    taille = 0
    indsend = 0

    if SIZE > 1:
        for i in range(len(halos.neigh[0])):
            scount[halos.neigh[0][i]] = halos.neigh[1][i]

        for i in range(SIZE):
            if i > 0:
                sdepl[i] = sdepl[i-1] + scount[i-1]

        rcount = COMM.alltoall(scount)

        for i in range(SIZE):
            if i > 0:
                rdepl[i] = rdepl[i-1] + rcount[i-1]

        for i in range(SIZE):
            taille += rcount[i]

        taille = int(taille)

        indsend = np.zeros(0, dtype=int)
        for i in range(len(halos.halosint)):
            indsend = np.append(indsend, cells.globtoloc[halos.halosint[i]])
            
    else:
        indsend = np.zeros(1, dtype=int)

    rcount = np.asarray(rcount)

    return scount, sdepl, rcount, rdepl, taille, indsend


def update_haloghost_info_2d(nodes, cells, halos ):
    
    nbnodes = len(nodes.name)
    
    ghostcenter = {}
    
    scount_node = np.zeros(SIZE, dtype=int)
    sdepl_node = np.zeros(SIZE, dtype=int)
    
    rcount_node = np.zeros(SIZE, dtype=int)
    rdepl_node = np.zeros(SIZE, dtype=int)
    taille = 0
    
    taille_node_ghost = np.zeros(SIZE, dtype=int)
    count = 6
    
    import collections
    if SIZE > 1:
        for i in range(nbnodes):
            if nodes.name[i] == 10:
                for j in range(len(nodes.ghostcenter[i])):
                    # if nodes.ghostcenter[i][j][-1] != -1:
                    for k in nodes.nparts[nodes.loctoglob[i]]:
                        if k != RANK:
                            ghostcenter.setdefault(k, []).append([
                                nodes.loctoglob[i], nodes.ghostcenter[i][j][0],
                                nodes.ghostcenter[i][j][1],
                                cells.loctoglob[nodes.ghostcenter[i][j][2]],
                                nodes.ghostcenter[i][j][3], nodes.ghostcenter[i][j][4]])

                            taille_node_ghost[k] += 1
        
        # print(ghostcenter)
        # import sys; sys.exit()
        ghostcenter = collections.OrderedDict(sorted(ghostcenter.items()))
        
        for i in range(len(halos.neigh[0])):
            scount_node[halos.neigh[0][i]] = taille_node_ghost[halos.neigh[0][i]]
        
        for i in range(1, SIZE):
            sdepl_node[i] = sdepl_node[i-1] + scount_node[i-1]
        
        rcount_node = COMM.alltoall(scount_node)
        
        for i in range(1, SIZE):
            rdepl_node[i] = rdepl_node[i-1] + rcount_node[i-1]
        
        for i in range(SIZE):
            taille += rcount_node[i]
    
        sendbuf = []
        
        for i, j in ghostcenter.items():
            sendbuf.extend(j)
        
        sendbuf = np.asarray(sendbuf)
        
        ghostcenter_halo = np.ones((taille, count))
        
        
        type_ligne = MPI.DOUBLE.Create_contiguous(count)
        type_ligne.Commit()
        
        s_msg = r_msg = 0
        s_msg = [sendbuf, (scount_node, sdepl_node), type_ligne]
        r_msg = [ghostcenter_halo, (rcount_node, rdepl_node), type_ligne]
        
        COMM.Alltoallv(s_msg, r_msg)
        
        recvbuf = {}
        for i in range(len(ghostcenter_halo)):
            recvbuf.setdefault(ghostcenter_halo[i][0], []).append([ghostcenter_halo[i][1], ghostcenter_halo[i][2], 
                                                                   ghostcenter_halo[i][3], ghostcenter_halo[i][4],
                                                                   ghostcenter_halo[i][5]])
        for i in range(nbnodes):
            if nodes.name[i] == 10:
                if recvbuf.get(nodes.loctoglob[i]):
                    nodes.haloghostcenter[i].extend(recvbuf[nodes.loctoglob[i]])
    
    maxGhostCell = 0
    for i in range(len(nodes.name)):
        maxGhostCell = max(maxGhostCell, len(nodes.ghostcenter[i]))
    
    for i in range(len(nodes.name)):
        iterator = maxGhostCell - len(nodes.ghostcenter[i])
        for k in range(iterator):
            nodes.ghostcenter[i].append([-1., -1., -1. , -1., -1])
    
        if len(nodes.ghostcenter[i]) == 0 :
                nodes.ghostcenter[i].append([-1, -1., -1., -1., -1])
                
    for i in range(len(nodes.name)):
        iterator = maxGhostCell - len(nodes.haloghostcenter[i])
        for k in range(iterator):
            nodes.haloghostcenter[i].append([-1., -1., -1., -1, -1])
        
        if len(nodes.haloghostcenter[i]) == 0 :
            nodes.haloghostcenter[i].append([-1, -1., -1., -1, -1])   
            
    #local halo index of haloext
    haloexttoind = {}
    for i in range(nbnodes):
        if nodes.name[i] == 10:
            for j in range(nodes.halonid[i][-1]):
                haloexttoind[halos.halosext[nodes.halonid[i][j]][0]] = nodes.halonid[i][j]
    # maxsize = 0
    cmpt = 0
    for i in range(nbnodes):
        if nodes.name[i] == 10:
            for j in range(len(nodes.haloghostcenter[i])):
                if nodes.haloghostcenter[i][j][-1] != -1:
                    nodes.haloghostcenter[i][j][-3] = haloexttoind[int(nodes.haloghostcenter[i][j][-3])]
                    # maxsize = int(max(maxsize, nodes.haloghostcenter[i][j][-1]+1))
                    nodes.haloghostcenter[i][j][-1] = cmpt
                    cmpt = cmpt + 1
            
            
    
    maxsize = cmpt
        
    nodes.ghostcenter     = np.asarray(nodes.ghostcenter)
    nodes.haloghostcenter = np.asarray(nodes.haloghostcenter)
        
    return maxsize

def update_haloghost_info_3d(nodes, cells, halos):
    
    nbnodes = len(nodes.name)
        
    ghostcenter = {}
    
    scount_node = np.zeros(SIZE, dtype=int)
    sdepl_node = np.zeros(SIZE, dtype=int)
    
    rcount_node = np.zeros(SIZE, dtype=int)
    rdepl_node = np.zeros(SIZE, dtype=int)
    taille = 0
    
    taille_node_ghost = np.zeros(SIZE, dtype=int)
    count = 7
    
    import collections
    maxsize = 0
    if SIZE > 1:
        for i in range(nbnodes):
            if nodes.name[i] == 10:
                for j in range(len(nodes.ghostcenter[i])):
                    for k in nodes.nparts[nodes.loctoglob[i]]:
                        if k != RANK:
                            ghostcenter.setdefault(k, []).append([
                                nodes.loctoglob[i], nodes.ghostcenter[i][j][0],
                                nodes.ghostcenter[i][j][1], nodes.ghostcenter[i][j][2], 
                                cells.loctoglob[nodes.ghostcenter[i][j][3]],
                                nodes.ghostcenter[i][j][4], nodes.ghostcenter[i][j][5]])
                            taille_node_ghost[k] += 1
        
        #print(maxsize)
        ghostcenter = collections.OrderedDict(sorted(ghostcenter.items()))
        
        for i in range(len(halos.neigh[0])):
            scount_node[halos.neigh[0][i]] = taille_node_ghost[halos.neigh[0][i]]
        
        for i in range(1, SIZE):
            sdepl_node[i] = sdepl_node[i-1] + scount_node[i-1]
        
        rcount_node = COMM.alltoall(scount_node)
        
        for i in range(1, SIZE):
            rdepl_node[i] = rdepl_node[i-1] + rcount_node[i-1]
        
        for i in range(SIZE):
            taille += rcount_node[i]
    
        sendbuf = []
        
        for i, j in ghostcenter.items():
            sendbuf.extend(j)
            
        sendbuf = np.asarray(sendbuf)
        
        ghostcenter_halo = np.ones((taille, count))
        
        
        type_ligne = MPI.DOUBLE.Create_contiguous(count)
        type_ligne.Commit()
        
        s_msg = r_msg = 0
        s_msg = [sendbuf, (scount_node, sdepl_node), type_ligne]
        r_msg = [ghostcenter_halo, (rcount_node, rdepl_node), type_ligne]
        
        COMM.Alltoallv(s_msg, r_msg)
        
        recvbuf = {}
        for i in range(len(ghostcenter_halo)):
            recvbuf.setdefault(ghostcenter_halo[i][0], []).append([ghostcenter_halo[i][1], ghostcenter_halo[i][2], 
                                                                   ghostcenter_halo[i][3], ghostcenter_halo[i][4],
                                                                   ghostcenter_halo[i][5], ghostcenter_halo[i][6]])
        for i in range(nbnodes):
            if nodes.name[i] == 10:
                if recvbuf.get(nodes.loctoglob[i]):
                    nodes.haloghostcenter[i].extend(recvbuf[nodes.loctoglob[i]])
    
    maxGhostCell = 0
    for i in range(len(nodes.name)):
        maxGhostCell = max(maxGhostCell, len(nodes.ghostcenter[i]))
    
    for i in range(len(nodes.name)):
        iterator = maxGhostCell - len(nodes.ghostcenter[i])
        for k in range(iterator):
            nodes.ghostcenter[i].append([-1., -1., -1., -1., -1., -1.])
    
        if len(nodes.ghostcenter[i]) == 0 :
                nodes.ghostcenter[i].append([-1, -1., -1.,-1., -1., -1.])
                
    maxGhostCell = 0
    for i in range(len(nodes.name)):
        maxGhostCell = max(maxGhostCell, len(nodes.haloghostcenter[i]))
        
    for i in range(len(nodes.name)):
        iterator = maxGhostCell - len(nodes.haloghostcenter[i])
        for k in range(iterator):
            nodes.haloghostcenter[i].append([-1., -1., -1., -1., -1., -1.])
    
        if len(nodes.haloghostcenter[i]) == 0 :
                nodes.haloghostcenter[i].append([-1, -1., -1.,-1., -1., -1.])
                
    #local halo index of haloext
    haloexttoind = {}
    for i in range(nbnodes):
        if nodes.name[i] == 10:
            for j in range(nodes.halonid[i][-1]):
                haloexttoind[halos.halosext[nodes.halonid[i][j]][0]] = nodes.halonid[i][j]
    cmpt = 0
    for i in range(nbnodes):
        if nodes.name[i] == 10:
            for j in range(len(nodes.haloghostcenter[i])):
                if nodes.haloghostcenter[i][j][-1] != -1:
                    nodes.haloghostcenter[i][j][-3] = haloexttoind[int(nodes.haloghostcenter[i][j][-3])]
                    nodes.haloghostcenter[i][j][-1] = cmpt
                    cmpt = cmpt + 1
    maxsize = cmpt#int(max(maxsize, nodes.haloghostcenter[i][j][-1]+1))
    
    nodes.ghostcenter     = np.asarray(nodes.ghostcenter)
    nodes.haloghostcenter = np.asarray(nodes.haloghostcenter)
    
    return maxsize
