from pyccel.decorators import stack_array, pure


@pure
@stack_array('norm', 'snom')
def create_info_2dfaces(cellid:'int[:,:]', nodeid:'int[:,:]', namen:'int[:]', vertex:'double[:,:]', 
                        centerc:'double[:,:]', nbfaces:'int', normalf:'double[:,:]', mesuref:'double[:]',
                        centerf:'double[:,:]', namef:'int[:]'):
    
    from numpy import double, zeros, sqrt
    
    norm   = zeros(3, dtype=double)
    snorm  = zeros(3, dtype=double)
    
    
    #Faces aux bords (1,2,3,4), Faces à l'interieur 0    A VOIR !!!!!
    for i in range(nbfaces):
        if (cellid[i][1] == -1 and cellid[i][1] != -10):
            if namen[nodeid[i][0]] == namen[nodeid[i][1]]:
                namef[i] = namen[nodeid[i][0]]
            elif ((namen[nodeid[i][0]] == 3 and namen[nodeid[i][1]] != 0) or
                    (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] == 3)):
                namef[i] = 3
            elif ((namen[nodeid[i][0]] == 4 and namen[nodeid[i][1]] != 0) or
                    (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] == 4)):
                namef[i] = 4
            elif ((namen[nodeid[i][0]] == 7 and namen[nodeid[i][1]] != 0) or
                    (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] == 7)):
                namef[i] = 7
            elif ((namen[nodeid[i][0]] == 8 and namen[nodeid[i][1]] != 0) or
                    (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] == 8)):
                namef[i] = 8    
            
            else:
                namef[i] = 100
        
        norm[0] = vertex[nodeid[i][0]][1] - vertex[nodeid[i][1]][1]
        norm[1] = vertex[nodeid[i][1]][0] - vertex[nodeid[i][0]][0]
    
        centerf[i][:] = 0.5 * (vertex[nodeid[i][0]][0:3] + vertex[nodeid[i][1]][0:3])
    
        snorm[:] = centerc[cellid[i][0]][:] - centerf[i][:]
    
        if (snorm[0] * norm[0] + snorm[1] * norm[1]) > 0:
            normalf[i][:] = -1*norm[:]
        else:
            normalf[i][:] = norm[:]

        mesuref[i] = sqrt(normalf[i][0]**2 + normalf[i][1]**2)

@pure   
@stack_array('norm', 'snom', 'u', 'v')
def create_info_3dfaces(cellid:'int[:,:]', nodeid:'int[:,:]', namen:'int[:]', vertex:'double[:,:]', 
                        centerc:'double[:,:]', nbfaces:'int', normalf:'double[:,:]', mesuref:'double[:]',
                        centerf:'double[:,:]', namef:'int[:]'):
    
    from numpy import double, zeros, sqrt
    
    norm   = zeros(3, dtype=double)
    snorm  = zeros(3, dtype=double)
    u      = zeros(3, dtype=double)
    v      = zeros(3, dtype=double)
    
    for i in range(nbfaces):
        if (cellid[i][1] == -1 ):
            if namen[nodeid[i][0]] == namen[nodeid[i][1]] and namen[nodeid[i][0]] == namen[nodeid[i][2]] :
                namef[i] = namen[nodeid[i][0]]
            
            elif ((namen[nodeid[i][0]] == 3 and namen[nodeid[i][1]] != 0 and namen[nodeid[i][2]] != 0) or
                  (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] == 3 and namen[nodeid[i][2]] != 0) or 
                  (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] != 0 and namen[nodeid[i][2]] != 3)):
                    namef[i] = 3
            
            elif ((namen[nodeid[i][0]] == 4 and namen[nodeid[i][1]] != 0 and namen[nodeid[i][2]] != 0) or
                (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] == 4 and namen[nodeid[i][2]] != 0) or 
                (namen[nodeid[i][0]] != 0 and namen[nodeid[i][1]] != 0 and namen[nodeid[i][2]] != 4)):
                namef[i] = 4
            else:
                namef[i] = 100
                
        u[:] = vertex[nodeid[i][1]][0:3]-vertex[nodeid[i][0]][0:3]
        v[:] = vertex[nodeid[i][2]][0:3]-vertex[nodeid[i][0]][0:3]
        
        norm[0] = 0.5*(u[1]*v[2] - u[2]*v[1])
        norm[1] = 0.5*(u[2]*v[0] - u[0]*v[2])
        norm[2] = 0.5*(u[0]*v[1] - u[1]*v[0])
    
        centerf[i][:] = 1./3 * (vertex[nodeid[i][0]][:3] + vertex[nodeid[i][1]][:3] + vertex[nodeid[i][2]][:3])
    
        snorm[:] = centerc[cellid[i][0]][:] - centerf[i][:]
    
        if (snorm[0] * norm[0] + snorm[1] * norm[1] + snorm[2] * norm[2]) > 0:
            normalf[i][:] = -1*norm[:]
        else:
            normalf[i][:] = norm[:]
    
        mesuref[i] = sqrt(normalf[i][0]**2 + normalf[i][1]**2 + normalf[i][2]**2)
        
@pure
def Compute_2dcentervolumeOfCell(nodeid:'int[:,:]', vertex:'double[:,:]', nbelements:'int',
                                 center:'double[:,:]', volume:'double[:]'):
    
    #calcul du barycentre et volume
    for i in range(nbelements):
        s_1 = nodeid[i][0]
        s_2 = nodeid[i][1]
        s_3 = nodeid[i][2]

        x_1 = vertex[s_1][0]; y_1 = vertex[s_1][1]
        x_2 = vertex[s_2][0]; y_2 = vertex[s_2][1]
        x_3 = vertex[s_3][0]; y_3 = vertex[s_3][1]

        center[i][0] = 1./3 * (x_1 + x_2 + x_3); center[i][1] = 1./3*(y_1 + y_2 + y_3); center[i][2] =  0.
        volume[i] = (1./2) * abs((x_1-x_2)*(y_1-y_3)-(x_1-x_3)*(y_1-y_2))
      
        var1 = (x_2-x_1)*(y_3-y_1)-(y_2-y_1)*(x_3-x_1)
        if var1 < 0:
            nodeid[i][0] = s_1;   nodeid[i][1] = s_3; nodeid[i][2] = s_2

@pure
@stack_array('u', 'v','w')
def Compute_3dcentervolumeOfCell(nodeid:'int[:,:]', vertex:'double[:,:]', nbelements:'int',
                                 center:'double[:,:]', volume:'double[:]'):
    
    from numpy import zeros, fabs
    wedge = zeros(3)
    u = zeros(3)
    v = zeros(3)
    w = zeros(3)
    
    #calcul du barycentre et volume
    for i in range(nbelements):
        
        s_1 = nodeid[i][0]
        s_2 = nodeid[i][1]
        s_3 = nodeid[i][2]
        s_4 = nodeid[i][3]
        
        x_1 = vertex[s_1][0]; y_1 = vertex[s_1][1]; z_1 = vertex[s_1][2]
        x_2 = vertex[s_2][0]; y_2 = vertex[s_2][1]; z_2 = vertex[s_2][2]
        x_3 = vertex[s_3][0]; y_3 = vertex[s_3][1]; z_3 = vertex[s_3][2]
        x_4 = vertex[s_4][0]; y_4 = vertex[s_4][1]; z_4 = vertex[s_4][2]
        
        center[i][0] = 1./4*(x_1 + x_2 + x_3 + x_4) 
        center[i][1] = 1./4*(y_1 + y_2 + y_3 + y_4)
        center[i][2] = 1./4*(z_1 + z_2 + z_3 + z_4)
        
        u[:] = vertex[s_2][:]-vertex[s_1][:]
        v[:] = vertex[s_3][:]-vertex[s_1][:]
        w[:] = vertex[s_4][:]-vertex[s_1][:]
        
        wedge[0] = v[1]*w[2] - v[2]*w[1]
        wedge[1] = v[2]*w[0] - v[0]*w[2]
        wedge[2] = v[0]*w[1] - v[1]*w[0]
        
        volume[i] = 1./6*fabs(u[0]*wedge[0] + u[1]*wedge[1] + u[2]*wedge[2]) 
        
@pure
def create_cellsOfFace(faceid:'int[:,:]', nbelements:'int', nbfaces:'int', cellid:'int[:,:]', dim:'int'):

    for i in range(nbelements):
        for j in range(dim+1):
            if cellid[faceid[i][j]][0] == -1 :
                cellid[faceid[i][j]][0] = i

            if cellid[faceid[i][j]][0] != i:
                cellid[faceid[i][j]][0] = cellid[faceid[i][j]][0]
                cellid[faceid[i][j]][1] = i

@pure
def create_2dfaces(nodeidc:'int[:,:]', nbelements:'int', faces:'int[:,:]',
                   cellf:'int[:,:]'):
   
    #Create 2d faces
    k = 0
    for i in range(nbelements):
        faces[k][0]   = nodeidc[i][0]; faces[k][1]   = nodeidc[i][1]
        faces[k+1][0] = nodeidc[i][1]; faces[k+1][1] = nodeidc[i][2]
        faces[k+2][0] = nodeidc[i][2]; faces[k+2][1] = nodeidc[i][0]
        cellf[i][0] = k; cellf[i][1] = k+1; cellf[i][2] = k+2
        k = k+3

@pure
def create_cell_faceid(nbelements:'int', oldTonewIndex:'int[:]', cellf:'int[:,:]', 
                       faceid:'int[:,:]', dim:'int'):
    
    for i in range(nbelements):
        for j in range(dim+1):
            faceid[i][j] = oldTonewIndex[cellf[i][j]]
        
@pure
def create_3dfaces(nodeidc:'int[:,:]', nbelements:'int',faces:'int[:,:]',
                   cellf:'int[:,:]'):
    #Create 3d faces
    k = 0
    
    for i in range(nbelements):
        faces[k][0]   = nodeidc[i][0]; faces[k][1]   = nodeidc[i][1]; faces[k][2]   = nodeidc[i][2]
        faces[k+1][0] = nodeidc[i][2]; faces[k+1][1] = nodeidc[i][3]; faces[k+1][2] = nodeidc[i][0]
        faces[k+2][0] = nodeidc[i][0]; faces[k+2][1] = nodeidc[i][1]; faces[k+2][2] = nodeidc[i][3]
        faces[k+3][0] = nodeidc[i][3]; faces[k+3][1] = nodeidc[i][1]; faces[k+3][2] = nodeidc[i][2]
        cellf[i][0]  = k; cellf[i][1] = k+1; cellf[i][2] = k+2; cellf[i][3] = k+3
        k = k+4
    
@pure
@stack_array('ss', 'G', 'c')
def create_NormalFacesOfCell(centerc:'double[:,:]', centerf:'double[:,:]', faceid:'int[:,:]', normal:'double[:,:]',
                             nbelements:'int', nf:'double[:,:,:]', dim:'int'):
    
    from numpy import zeros
    ss = zeros(3)
    G  = zeros(3)
    c  = zeros(3)
    
    #compute the outgoing normal faces for each cell
    for i in range(nbelements):
        G[:] = centerc[i][:]

        for j in range(dim+1):
            f = faceid[i][j]
            c[:] = centerf[f][:]

            if ((G[0]-c[0])*normal[f][0] + (G[1]-c[1])*normal[f][1] + (G[2]-c[2])*normal[f][2]) < 0.:
                ss[:] = normal[f][:]
            else:
                ss[:] = -1.0*normal[f][:]
                
            nf[i][j][:] = ss[:]
