#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

if(argc != 5){
    printf("Usage %s, n_iters, n, m, k\n", argv[0]);
    return -1;
}

int n_iters = atoi(argv[1]);
int n = atoi(argv[2]);
int m = atoi(argv[3]);
int k = atoi(argv[4]);

if(!(m>n)){
    printf("m must be strictly greater than n\n");
    return 69;
}

int rank, num_proc, left, right;
MPI_Init(&argc, &argv);
MPI_Status status;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

if(rank==0)
if(m != (k*num_proc)){
    printf("m != k*num_proc");
    return 1;
}

MPI_Comm CIRC_COMM;
int dims[1];
int periods[1];
int reorder = 0;
dims[0] = num_proc;
periods[0] = 1;

MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &CIRC_COMM);
MPI_Cart_shift(CIRC_COMM, 0, -1, &right, &left);


MPI_Datatype rtype;
MPI_Type_contiguous(n*m, MPI_INT, &rtype);
MPI_Type_commit(&rtype);

MPI_Datatype ctype;
MPI_Type_vector(k, 1, n, MPI_INT, &ctype);
MPI_Type_commit(&ctype);


int *localA = (int *)malloc(k*n*sizeof(int));
int *At = (int *)malloc(n*k*sizeof(int));
int *localb = (int *)malloc(n*sizeof(int));
int *localw = (int *)malloc(k*sizeof(int));
for(int i = 0; i < k*n; i++){
    localA[i] = i;//rand() % n;
}
for(int i = 0; i < n; i++){
    localb[i] = rand() % k;
}

// for(int i = 0; i < k; i++){
//     printf("Rank %d[", rank);
//     for(int j = 0; j < n; j++){
//         printf("%d ", localA[i*n + j]);
//     }
//     printf("]\n");
// }
// Transpose matrix

for(int i = 0; i < k; i++){
    for(int j = i + 1; j < n; j++){
        At[i*n + j] = localA[j*n + i];
        At[j*n + i] = localA[i*n + j];;
    }
}


// for(int i = 0; i < k; i++){
//     printf("Rank %d[", rank);
//     for(int j = 0; j < n; j++){
//         printf("%d ", localA[i*n + j]);
//     }
//     printf("]\n");
// }

MPI_Finalize();

return 0;
}