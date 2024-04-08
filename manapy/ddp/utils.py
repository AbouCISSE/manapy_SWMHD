#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 21:15:29 2021

@author: kissami
"""
from numpy import zeros, asarray, double, int64, dot, cross
from numba import njit

@njit
def oriente_3dfacenodeid(nodeid:'int[:,:]', normal:'double[:,:]', vertex:'double[:,:]'):
    
    nbfaces = len(nodeid)
    v1 = zeros(3)
    v2 = zeros(3)
    
    for i in range(nbfaces):
        n1 = nodeid[i][0]; n2 = nodeid[i][1]; n3 = nodeid[i][2]
        s1 = vertex[n1][0:3]; s2 = vertex[n2][0:3]; s3 = vertex[n3][0:3];
        
        v1[:] = s2[:] - s1[:] 
        v2[:] = s3[:] - s1[:]
        
        if dot(cross(v1, v2), normal[i]) < 0:
            nodeid[i][1] = n3; nodeid[i][2] = n2
            
    return nodeid

@njit 
def create_NeighborCellByFace(faceid:'int[:,:]', cellid:'int[:,:]', nbelements:'int', dim:'int'):
    
    cellfid = [[i for i in range(0)] for i in range(nbelements)]
    #Création des 3/4 triangles voisins par face
    for i in range(nbelements):
        for j in range(dim+1):
            f = faceid[i][j]
            if cellid[f][1] != -1:
                if i == cellid[f][0]:
                    cellfid[i].append(cellid[f][1])
                else:
                    cellfid[i].append(cellid[f][0])
            else:
                cellfid[i].append(-1)
    cellfid = asarray(cellfid, dtype=int64)
    
    return cellfid

@njit 
def create_2doppNodeOfFaces(nodeidc:'int[:,:]', faceidc:'int[:,:]', nodeidf:'int[:,:]', nbelements:'int', nbfaces:'int'):
    #    #TODO improve the opp node creation
    oppnodeid = [[i for i in range(0)] for i in range(nbfaces)]
   
    for i in range(nbelements):
        f1 = faceidc[i][0]; f2 = faceidc[i][1]; f3 = faceidc[i][2]
        n1 = nodeidc[i][0]; n2 = nodeidc[i][1]; n3 = nodeidc[i][2] 
       
        if n1 not in nodeidf[f1] :
            oppnodeid[f1].append(n1)
        if n1 not in nodeidf[f2] :
            oppnodeid[f2].append(n1)
        if n1 not in nodeidf[f3] :
            oppnodeid[f3].append(n1)
        
        if n2 not in nodeidf[f1] :
            oppnodeid[f1].append(n2)
        if n2 not in nodeidf[f2] :
            oppnodeid[f2].append(n2)
        if n2 not in nodeidf[f3] :
            oppnodeid[f3].append(n2)
        
        if n3 not in nodeidf[f1] :
            oppnodeid[f1].append(n3)
        if n3 not in nodeidf[f2] :
            oppnodeid[f2].append(n3)
        if n3 not in nodeidf[f3] :
            oppnodeid[f3].append(n3)
        
    for i in range(nbfaces):
        if len(oppnodeid[i]) < 2:
            oppnodeid[i].append(-1)
    
    oppnodeid = asarray(oppnodeid, dtype=int64)      
            
    return oppnodeid

@njit 
def create_3doppNodeOfFaces(nodeidc:'int[:,:]', faceidc:'int[:,:]', nodeidf:'int[:,:]', nbelements:'int', nbfaces:'int'):
    
    #TODO improve the opp node creation
    oppnodeid = [[i for i in range(0)] for i in range(nbfaces)]
    for i in range(nbelements):
        f1 = faceidc[i][0]; f2 = faceidc[i][1]; f3 = faceidc[i][2]; f4 = faceidc[i][3]; 
        n1 = nodeidc[i][0]; n2 = nodeidc[i][1]; n3 = nodeidc[i][2]; n4 = nodeidc[i][3]; 
        
        if n1 not in nodeidf[f1] :
            oppnodeid[f1].append(n1)
        if n1 not in nodeidf[f2] :
            oppnodeid[f2].append(n1)
        if n1 not in nodeidf[f3] :
            oppnodeid[f3].append(n1)
        if n1 not in nodeidf[f4] :
            oppnodeid[f4].append(n1)
        
        if n2 not in nodeidf[f1] :
            oppnodeid[f1].append(n2)
        if n2 not in nodeidf[f2] :
            oppnodeid[f2].append(n2)
        if n2 not in nodeidf[f3] :
            oppnodeid[f3].append(n2)
        if n2 not in nodeidf[f4] :
            oppnodeid[f4].append(n2)
        
        if n3 not in nodeidf[f1] :
            oppnodeid[f1].append(n3)
        if n3 not in nodeidf[f2] :
            oppnodeid[f2].append(n3)
        if n3 not in nodeidf[f3] :
            oppnodeid[f3].append(n3)
        if n3 not in nodeidf[f4] :
            oppnodeid[f4].append(n3)
        
        if n4 not in nodeidf[f1] :
            oppnodeid[f1].append(n4)
        if n4 not in nodeidf[f2] :
            oppnodeid[f2].append(n4)
        if n4 not in nodeidf[f3] :
            oppnodeid[f3].append(n4)
        if n4 not in nodeidf[f4] :
            oppnodeid[f4].append(n4)
        
            
    for i in range(nbfaces):
        if len(oppnodeid[i]) < 2:
            oppnodeid[i].append(-1)
            
    oppnodeid = asarray(oppnodeid, dtype=int64)   
    
    return oppnodeid


@njit 
def create_node_cellid(nodeid:'int[:,:]', vertex:'double[:,:]', nbelements:'int', nbnodes:'int', dim:'int'):
    
    tmp = [[i for i in range(0)] for i in range(nbnodes)]
    longn = zeros(nbnodes, dtype=int64)
    
    for i in range(nbelements):
        for j in range(dim+1):
            tmp[nodeid[i][j]].append(i)
            longn[nodeid[i][j]] = longn[nodeid[i][j]] + 1
    
    longc = zeros(nbelements, dtype=int64)
    tmp2 = [[i for i in range(0)] for i in range(nbelements)]
    
    for i in range(nbelements):
        for j in range(dim+1):
            for k in range(len(tmp[nodeid[i][j]])):
                if (tmp[nodeid[i][j]][k] not in tmp2[i] and  tmp[nodeid[i][j]][k] != i):
                    tmp2[i].append(tmp[nodeid[i][j]][k])
                    longc[i] = longc[i] + 1

    maxlongn = int(max(longn))
    cellid = [[-1 for i in range(maxlongn)] for i in range(nbnodes)]
    
    for i in range(len(tmp)):
        for j in range(len(tmp[i])):
            cellid[i][j] = tmp[i][j]
            
    for i in range(nbnodes):
        cellid[i].append(longn[i])
    cellid = asarray(cellid, dtype=int64)
    
    maxlongc = int(max(longc))
    cellnid = [[-1 for i in range(maxlongc)] for i in range(nbelements)]
    
    for i in range(len(tmp2)):
        for j in range(len(tmp2[i])):
            cellnid[i][j] = tmp2[i][j]
    
    for i in range(nbelements):
        cellnid[i].append(longc[i])
   
    cellnid = asarray(cellnid, dtype=int64)   
    
    
    return cellid, cellnid
