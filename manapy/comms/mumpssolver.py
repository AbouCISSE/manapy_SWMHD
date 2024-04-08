"""
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

#mpicc -fPIC -shared -o mumps.so  mumps.c /usr/lib/x86_64-linux-gnu/libmumps_common-5.1.2.so /usr/lib/x86_64-linux-gnu/libdmumps-5.1.2.so
try:
    COMMFILE = os.environ['COMMFILE']
except:
    BASE_DIR = os.path.dirname(os.path.realpath(__file__))
    COMMFILE = os.path.join(BASE_DIR, '.')
filename = os.path.join(COMMFILE, "mumps.so")

commFunc = CDLL(filename)

commFunc.mumps_solve.argtypes = [POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_double), 
                                 c_int, c_int, c_int] 


def mumps_solve(a_loc, irn_loc, jcn_loc, rhs, nz_loc, n, rank):
    
    commFunc.mumps_solve(a_loc.ctypes.data_as(POINTER(c_double)), irn_loc.ctypes.data_as(POINTER(c_int)), 
                        jcn_loc.ctypes.data_as(POINTER(c_int)), rhs.ctypes.data_as(POINTER(c_double)), 
                        nz_loc, n, rank)
    
