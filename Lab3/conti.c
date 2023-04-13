#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

int rank, num_proc, left, right;
MPI_Init(&argc, &argv);
MPI_Status status;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

MPI_Comm CIRC_COMM;
int dims[1];
int periods[1];
int reorder = 0;
dims[0] = num_proc;
periods[0] = 1;

MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &CIRC_COMM);
MPI_Cart_shift(CIRC_COMM, 0, -1, &right, &left);

int n = 10;
int k = n/num_proc;
int m = 2;

MPI_Datatype rtype;
MPI_Type_contiguous(n*m, MPI_INT, &rtype);
MPI_Type_commit(&rtype);

MPI_Datatype ctype;
MPI_Type_vector(k, 1, n, MPI_INT, &ctype);
MPI_Type_commit(&ctype);


int *globarray;
int *localarray = (int *)malloc(k*n*sizeof(int));
if(rank == 0){
    globarray = (int *)malloc(n*n*sizeof(int));
  for(int i = 0; i < n*n; i++){
  globarray[i] = i;
  }
  for(int j = 0; j<n;j++){
    for(int i = 0; i<n; i++){
        printf("%d ", globarray[i+j*n]);
    }
    printf("\n");
  }
}

MPI_Scatter(globarray, n*k, MPI_INT, localarray, n*k, MPI_INT, 0, CIRC_COMM);

// for(int i = 0; i<k; i++){
//     printf("Rank %d [", rank);
//     for(int j = 0; j < n; j++){
//         printf("%d ", localarray[i*n + j]);
//     }
//     printf("\n");
// }

// Send m last rows
// MPI_Sendrecv(&localarray[(k - m)*n], m*n, MPI_INT, right, 1, &localarray[(k - m)*n], 1, rtype, left, 1, CIRC_COMM, &status);

// Send last column
MPI_Sendrecv(&localarray[n-1], 1, ctype, right, 1, localarray, 1, ctype, left, 1, CIRC_COMM, &status);

for(int i = 0; i<k; i++){
    printf("Rank %d [", rank);
    for(int j = 0; j < n; j++){
        printf("%d ", localarray[i*n + j]);
    }
    printf("\n");
}


MPI_Finalize();

return 0;
}