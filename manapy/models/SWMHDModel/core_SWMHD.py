#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 23:39:55 2022

@author: kissami
"""

from manapy.ast.core import Variable

#from manapy.models.SWMHDModel.tools import term_coriolis_SW
from manapy.models.SWMHDModel.tools import (update_SWMHD, time_step_SWMHD, explicitscheme_convective_SWMHD,
                                               term_source_srnh_SWMHD, explicitscheme_stabilizater, cpsi_global)

from mpi4py import MPI


class ShallowWaterMHDModel():
    
    variables=["h", "hu", "hv", "hB1", "hB2", "PSI", "Z"]
    varbs = {}
    def __init__(self, domain=None, terms=None, parameters=None, scheme=None, order=None, comm=None, *args, **kwargs):
        
        if domain is None:
            raise ValueError("domain must be given")
        if order is None:
            order = 1
        if comm is None:
            comm = MPI.COMM_WORLD
            
        if terms is None:
            terms = ['source', 'stabilizater', 'coriolis', "convective"]
        else:
            terms.extend(["convective"])

        self.domain = domain
        self.terms = terms
        self.parameters = parameters
        self.grav = 1.0
        self.GLM = 100
        self.order = order
        self.comm = comm
        self.Dz = 0.        
        
        if "stabilizater" in self.terms:
            if "vepsilon" not in self.parameters:
                raise ValueError("epsilon number must be given")
            else:
                self.vepsilon  = self.parameters["vepsilon"]
        else:  
            self.vepsilon  = 0.        
        
        if "coriolis" in self.terms:
            if "CoriolisForce" not in self.parameters:
                raise ValueError("Coriolis Force number must be given")
            else:
                self.fc  = self.parameters["Manning"]
        else:  
            self.fc  = 0.
        
        for var in self.variables:
            self.varbs[var] = Variable(domain=self.domain, terms=self.terms, parameters=self.parameters, name=var)
        
        self.h   = self.varbs['h']
        self.hu  = self.varbs['hu']
        self.hv  = self.varbs['hv']
        self.hB1 = self.varbs['hB1']
        self.hB2 = self.varbs['hB2']
        self.PSI = self.varbs['PSI']
        self.Z   = self.varbs['Z']
        
        
    def update_values(self):
        for var in self.varbs.values():
            var.update_values()





    def update_cpsi_global(self, cfl=None, comm=None):
        cpsi=cpsi_global(self.h.cell, self.hu.cell, self.hv.cell, self.hB1.cell, self.hB2.cell, cfl, self.domain.faces._normal, self.domain.faces._mesure, 
                     self.domain.cells._volume, self.domain.cells._faceid)
   
        self.cpsi = self.comm.allreduce(cpsi, MPI.MAX)
        
    def update_explicitscheme_convective(self, order=None):
        
        if order == 2:
            self.h.compute_cell_gradient()

        explicitscheme_convective_SWMHD(self.h.convective, self.hu.convective, self.hv.convective, self.hB1.convective, self.hB2.convective,
        			      self.PSI.convective, self.Z.convective, 
                                     self.h.cell, self.hu.cell, self.hv.cell, self.hB1.cell, self.hB2.cell, self.PSI.cell, self.Z.cell,
                                     self.h.ghost, self.hu.ghost, self.hv.ghost, self.hB1.ghost, self.hB2.ghost, self.PSI.ghost, self.Z.ghost, self.h.halo, 
                                     self.hu.halo, self.hv.halo, self.hB1.halo, self.hB2.halo, self.PSI.halo, self.Z.halo,
                                     self.h.gradcellx, self.h.gradcelly, self.h.gradhalocellx, self.h.gradhalocelly,
                                     self.h.psi, self.h.psihalo, 
                                     self.domain.cells._center, self.domain.faces._center, self.domain.halos._centvol, 
                                     self.domain.faces._ghostcenter, self.domain.faces._cellid, self.domain.faces._mesure, 
                                     self.domain.faces._normal, self.domain.faces._halofid, self.domain.innerfaces, self.domain.halofaces, 
                                     self.domain.boundaryfaces, self.domain.periodicboundaryfaces, self.domain.cells.shift, self.order, self.cpsi, self.hB1.cell, self.hB2.cell)



#    def update_explicit_stabilizater(self):
#    
#        #self.h.compute_face_gradient()
#        self.h.interpolate_celltonode()
#        self.hu.compute_face_gradient_uv(self.h.cell, self.h.ghost, self.h.halo, self.h.node)
#        self.hv.compute_face_gradient_uv(self.h.cell, self.h.ghost, self.h.halo, self.h.node)
#
#        explicitscheme_stabilizater(self.hu.gradfacex, self.hu.gradfacey, self.domain.faces.cellid, 
#                                   self.domain.faces.normal, self.domain.faces.name, self.vepsilon, self.hu.stabilizater, self.h.cell)
#
#        explicitscheme_stabilizater(self.hv.gradfacex, self.hv.gradfacey, self.domain.faces.cellid, 
#                                   self.domain.faces.normal, self.domain.faces.name, self.vepsilon, self.hv.stabilizater, self.h.cell)

    def update_term_source(self, order=1):
        
        term_source_srnh_SWMHD(self.h.source, self.hu.source, self.hv.source, self.hB1.source, self.hB2.source, self.PSI.source, self.Z.source,
                            self.h.cell, self.hu.cell, self.hv.cell, self.Z.cell,
                            self.h.ghost, self.hu.ghost, self.hv.ghost, self.Z.ghost,
                            self.h.halo, self.hu.halo, self.hv.halo, self.Z.halo,
                            self.h.gradcellx, self.h.gradcelly, self.h.psi, self.h.gradhalocellx, self.h.gradhalocelly, 
                            self.h.psihalo, 
                            self.domain.cells.nodeid, self.domain.cells.faceid, self.domain.cells.cellfid, self.domain.faces.cellid,
                            self.domain.cells.center, self.domain.cells.nf, 
                            self.domain.faces.name, self.domain.faces.center, self.domain.halos.centvol,
                            self.domain.nodes.vertex, self.domain.faces.halofid, self.order)
        
 
    
    
    
    
    
    
    
    
#    def update_term_coriolis(self):
#        term_coriolis_SW(self.hu.cell, self.hv.cell, self.hu.coriolis, self.hv.coriolis, self.fc)
        
    def update_time_step(self, cfl=None, comm=None):
        dt_c=time_step_SWMHD(self.h.cell, self.hu.cell, self.hv.cell, self.hB1.cell, self.hB2.cell, cfl, self.domain.faces._normal, self.domain.faces._mesure, 
                     self.domain.cells._volume, self.domain.cells._faceid)
        
        self.time_step = self.comm.allreduce(dt_c, MPI.MIN)


    
    def update_solution(self, order = None):
        
        if order is None:
            order = 1
 

        #update the convectiveidus using explicit scheme
        self.update_explicitscheme_convective(order=order)
    
        #update the source using explicit scheme
        self.update_term_source(order=order)
        
        #update coriolis forces
        #self.update_term_coriolis()
        
        #update stabilizater term       
        #self.update_explicit_stabilizater()
        

        update_SWMHD(self.h.cell, self.hu.cell, self.hv.cell, self.hB1.cell, self.hB2.cell, self.PSI.cell, self.Z.cell ,
                  self.h.convective, self.hu.convective, self.hv.convective, self.hB1.convective, self.hB2.convective, self.PSI.convective, self.Z.convective,
                  self. h.source, self.hu.source, self.hv.source, self.hB1.source, self.hB2.source, self.PSI.source, self.Z.source,
                  self.time_step, self.domain.cells._volume, self.GLM, self.cpsi)

    
