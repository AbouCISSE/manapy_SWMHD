#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "dmumps_c.h"   
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]


void mumps_solve(double* a_loc, int* irn_loc, int* jcn_loc, double* rhs, int nz_loc, int n, int rank) {


DMUMPS_STRUC_C id;

  printf("%d %d\n", rank, n);
  id.job=JOB_INIT; id.par=1; id.sym=1; id.comm_fortran=USE_COMM_WORLD;
  dmumps_c(&id);


  if(rank==0){
  id.n = n;
  }
  printf("%d %d\n", rank, n);
  
  id.ICNTL(7)  = 5;
  id.ICNTL(18) = 3;
  id.nz_loc    = nz_loc;
  id.a_loc     = a_loc;
  id.irn_loc   = irn_loc;
  id.jcn_loc   = jcn_loc;
  
//  for (int i = 0; i < nz_loc; i++)
//      printf("%f %d %d %d %d\n", id.a_loc[i], id.irn_loc[i], id.jcn_loc[i], id.n, id.nz_loc);
  //id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=5;
  
  /*
  id.job=2;
  dmumps_c(&id);
*/
  if(rank==0){
    id.rhs=&rhs[0];
  }

  id.job=6;
  dmumps_c(&id);

  /* Terminate instance. */
  id.job=JOB_END;
  dmumps_c(&id);


}
/*
void mumps_solve(int rank, double* b){
  #define ICNTL(I) icntl[(I)-1]

  if(rank==0){
    id.rhs=&b[0];
  }
  id.job=3;
  dmumps_c(&id);  
}
*/
