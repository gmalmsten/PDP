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
    MPI_Request request;
    if(rank == 0){
      MPI_Isend(&a, 1, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &request);
      MPI_Irecv(&b, 1, MPI_DOUBLE, previous, 0, MPI_COMM_WORLD, &request);
      printf("Processor %d waiting for data\n", rank);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
      printf("Processor %d received data\n", rank);
    }
    else{
      
      MPI_Request request;
      MPI_Irecv(&b, 1, MPI_DOUBLE, previous, 0, MPI_COMM_WORLD, &request);
      printf("Processor %d waiting for data\n", rank);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
      printf("Processor %d received data, sending to %d\n", rank, next);
      MPI_Isend(&a, 1, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &request);
    }
    // printf("Processor %d got %f from processor %d\n", rank, b, previous);
  

  MPI_Finalize(); 

  return 0;
}
