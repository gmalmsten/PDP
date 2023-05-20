#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include "prop.h"
#include <string.h>

#define ROWS 15
#define COLS 7
#define T 100
#define b 20

#define PRODUCE_OUTPUT

void print_d_vec(double *vector, int lim){
    /*Print a vector consisting of lim doubles*/
    printf("[");
    for(int i = 0; i < lim; i++){
        printf("%lf, ", vector[i]);
    }
    printf("]\n");
}

void print_i_vec(int *vector, int lim){
    /*Print a vector consisting of lim integers*/
    printf("[");
    for(int i = 0; i < lim - 1; i++){
        printf("%d, ", vector[i]);
    }
    printf("%d]\n", vector[lim - 1]);
}

int cmp (const void *num1, const void *num2) {
    /*Comparison function for qsort*/
   return ( *(int*)num1 > *(int*)num2 );
}

int main(int argc, char *argv[]){

    if(argc != 2){
        printf("Usage %s N\n", argv[0]);
        return -1;
    }

    // Arguments 
    const int N = atoi(argv[1]);
    int rank, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);   
    const int local_N = N/num_proc;
    
    // Assumption: N is divisible by num_proc
    if(N%num_proc){
        if(rank == 0){
            printf("ERROR N%%p != 0\n");
            return -2;
        }
    }
    
    
    // Transition matrix
    const int P[ROWS*COLS] = {1, 0, 0, 0, 0, 0, 0,
                            -1, 0, 0, 0, 0, 0, 0,
                            -1, 0, 1, 0, 0, 0, 0,
                            0, 1, 0, 0, 0, 0, 0,
                            0, -1, 0, 0, 0, 0, 0,
                            0, -1, 0, 1, 0, 0, 0,
                            0, 0, -1, 0, 0, 0, 0,
                            0, 0, -1, 0, 1, 0, 0,
                            0, 0, 0, -1, 0, 0, 0,
                            0, 0, 0, -1, 0, 1, 0,
                            0, 0, 0, 0, -1, 0, 0,
                            0, 0, 0, 0, -1, 0, 1,
                            0, 0, 0, 0, 0, -1, 0,
                            1, 0, 0, 0, 0, 0, -1,
                            0, 0, 0, 0, 0, 0, -1};

    // Seed
    // time_t seed = time(NULL);
    int seed = 1;
    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    srand(seed + rank);

    
    int results[local_N];

    // Start timer
    double start_time = MPI_Wtime();
    for(int epoch = 0; epoch < local_N; epoch++){
        // Initialize new simulation
        double t = 0;
        int x[COLS] = {900, 900, 30, 330, 50, 270, 20};

        char timing = 0;  // Current sub time to store, 0 - 25, 1 - 50, 2 - 75, 3 - 100
        double sub_start_time = MPI_Wtime();
        while(t<T){
            // Compute w 
            double w[ROWS];
            prop(x, w);

            // Compute a0
            double a0 = 0;
            for(int i = 0; i < ROWS; i++){
                a0+=w[i];
            }
            if(a0<0){
                printf("ERROR a0 < 0\n");
                return -3;
            }

            double tau = -log(((double)rand()/RAND_MAX))/a0;

            // Find r
            double sum_r_prev = 0;
            double sum_r = w[0];
            double lim = a0*((double)rand()/RAND_MAX);
            int r = 0;
            while(sum_r < lim && sum_r_prev <= lim){   //-----------------------------------------------------------//
                r++;
                sum_r += w[r];
                sum_r_prev = sum_r;
                if(r>ROWS){
                    printf("Error: r exceeds the bounds of the w array\n");
                    return -1;
                }
            }


            // Update x
            for(int i = 0; i < COLS; i++){
                x[i] += P[r*COLS + i];
            }

            // Step time
            t+=tau;
        }
        // Store sub time and result
        results[epoch] = x[0];
    }

    // Calculate local and global min and max
    int global_min, global_max, local_min, local_max;

    local_min = results[0];
    local_max = results[0];
    for(int i = 1; i < local_N; i++)
    {
        int tmp = results[i];
        if(tmp < local_min){
            local_min = tmp;
        }
        if(tmp > local_max){
            local_max = tmp;
        }
    }


    // Reduce and broadcast global min and max to all processes
    MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);


    // Calculate bin size and the local counts in each bin
    int bin_size = (global_max - global_min)/b;
    int bins[b] = {0};


    // Sort the local results and count the number of elements in each bin
    qsort(results, local_N, sizeof(int), cmp);

    int bin = 0;    // Current bin
    for(int i = 0; i < local_N; i++){
        if(results[i] > global_min + bin_size*(bin+1)){
            bin++;
        }
        if(results[i] > global_min + bin_size*b) // Let last bin contain all elements larger than global_min + 20*bin_size
        {
            bins[b-1]++;
        }
        else
        {
            bins[bin]++;
        }
    }

    // Sum all local results in root process (0)
    int global_bins[b];
    MPI_Reduce(bins, global_bins, b, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Stop timer
    double local_time = MPI_Wtime() - start_time;
    double global_time;
    MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    

    // Produce output
    if(rank == 0)
    {
        #ifdef PRODUCE_OUTPUT
        printf("Bin %d [%d %d]\n", 1, global_min, global_min + bin_size);
        for(int i = 1; i < b - 1; i++)
        {
            printf("Bin %d (%d %d]\n", i+1, global_min + bin_size*i, global_min + bin_size*(i + 1));
        }
        printf("Bin %d (%d %d]\n", 20, global_min + bin_size*19, global_max);
        print_i_vec(global_bins, b);
        #endif

        printf("%lf\n", global_time);
    }


    MPI_Finalize();
    return 0;
}