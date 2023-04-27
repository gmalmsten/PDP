#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void comm_splitR(int local_max, int local_min, int level, MPI_Comm communicator);

int main(int argc, char **argv) {

	// Parallel environment
	MPI_Init(&argc, & argv);
	int rank, num_proc;
	int N = atoi(argv[1]);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    int m  = N/num_proc;
    if(N%num_proc){
        if(rank == 0){
            printf("Error N != m*num_proc\n");
            return -1;
        }
    }
    
    int *list_to_sort;
    int *local_list = (int *)malloc(m*sizeof(int));

    if(rank == 0){
        list_to_sort = (int *)malloc(N*sizeof(int));
        for(int i = 0; i < N; i++){
            list_to_sort[i] = rand()%(N/2);
        }
    }

    MPI_Scatter(list_to_sort, m, MPI_INT, local_list, m, MPI_INT, 0, MPI_COMM_WORLD);


    int index = 0;
    while(1){

        if(index>=m){
            index = N;
            break;
        }
        if(local_list[index] == 0){
            printf("rank %d Found 0 at %d\n", rank, index);
            break;
        }
        index ++;
    }

    int my_index = rank*m + index;
    int global_index;

    MPI_Reduce(&my_index, &global_index, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    if(rank == 0){
        int i = 0;
        while(1){
            if(i >= N){
                printf("NO ZEROS IN LIST\n");
                break;
            }
            if(list_to_sort[i] == 0){
                break;
            }
            i++;
        }
        printf("First zero in list at index %d, check %d\n", global_index, i);
    }

    

    MPI_Finalize();
    return 0;
}