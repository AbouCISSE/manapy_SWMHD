#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>


void communication(double* sendbuf, double *recvbuf, int* displs, int* sendcounts , 
                    int* rdispls, int* recvcounts, void *ptr) {

	MPI_Comm graph_comm = *((MPI_Comm*)ptr);

	MPI_Neighbor_alltoallv(sendbuf, sendcounts, displs, MPI_DOUBLE, recvbuf, recvcounts, rdispls, MPI_DOUBLE, graph_comm);

}
