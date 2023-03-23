/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int rank, size;
  double a, b;
  MPI_Status status;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

  int previous = (rank - 1 + size)%size;
  int next = (rank + 1)%size;
  
  a = 100.0 + (double) rank;  /* Different a on different processors */

  /* Exchange variable a, notice the send-recv order */

    MPI_Send(&a, 1, MPI_DOUBLE, next, 69, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, previous, 69, MPI_COMM_WORLD, &status);
    // MPI_Sendrecv(&a, 1, MPI_DOUBLE, next, 69, &b, 1, MPI_DOUBLE, previous, 69, MPI_COMM_WORLD, &status);
    printf("Processor %d got %f from processor %d\n", rank, b, previous);
  

  MPI_Finalize(); 

  return 0;
}
