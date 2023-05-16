#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include "prop.h"

#define ROWS 15
#define COLS 7
#define T 100
#define b 20

void print_d_vec(double *vector, int lim){
    printf("[");
    for(int i = 0; i < lim; i++){
        printf("%f, ", vector[i]);
    }
    printf("]\n");
}

void print_i_vec(int *vector, int lim){
    printf("[");
    for(int i = 0; i < lim - 1; i++){
        printf("%d, ", vector[i]);
    }
    printf("%d]\n", vector[lim - 1]);
}

int cmp (const void *num1, const void *num2) {
   return ( *(int*)num1 > *(int*)num2 );
}

int main(int argc, char *argv[]){

    if(argc != 2){
        printf("Usage %s N\n", argv[0]);
        return -1;
    }

    // Arguments / 
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
    
    
    int X[COLS*local_N];
    int results[local_N];
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

    double start_time = MPI_Wtime();
    
    for(int epoch = 0; epoch < local_N; epoch++){
        srand((unsigned int)(rank*local_N + epoch));
        // srand((unsigned int) time(NULL));
        // printf("Rank %d: %d\n", rank, rank*local_N + i);
        double t = 0;
        int x[COLS] = {900, 900, 30, 330, 50, 270, 20};
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
                printf("ERROR a0 < 0\nw: [");
                for(int i = 0; i < ROWS; i++){
                    printf("%lf ", w[i]);
                }
                printf("]\nx: [");
                for(int i = 0; i < COLS; i++){
                    printf("%d ", x[i]);
                }
                printf("]\n");
                return -3;
            }

            // Generate two random numbers
            double u1 = (double)rand()/RAND_MAX;
            double u2 = (double)rand()/RAND_MAX;

            double tau = -log(u1)/a0;

            // Find r
            double sum_r_prev = 0;
            double sum_r = w[0];
            double lim = a0*u2;
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
                if(x[i]<0){
                    x[i] = 0;
                    printf("Warning: x[%d] < 0\n", i);
                }
            }

            // Step time
            t+=tau;
        }
        // for(int j = 0; j < n; j++){
        //     X[i*n + j] = x[j];
        // }
        results[epoch] = x[0];
    }

    // Calculate local and global min and max
    int global_min, global_max, local_min, local_max;

    local_min = results[0];
    local_max = results[0];
    for(int i = 1; i < local_N; i++){
        if(results[i] < local_min){
            local_min = results[i];
        }
        if(results[i] > local_max){
            local_max = results[i];
        }
    }


    MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);


    // Calculate bin size and the local counts in each bin
    int bin_size = (global_max - global_min)/b;
    int bins[b];
    for(int i = 0; i < b; i++){
        bins[i] = 0;
    }
    qsort(results, local_N, sizeof(int), cmp);

    int bin = 0;    // Current bin
    for(int i = 0; i < local_N; i++){
        if(results[i] > global_min + bin_size*(bin+1)){
            bin++;
        }
        bins[bin]++;
    }

    int global_bins[b];

    MPI_Reduce(bins, global_bins, b, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double local_time = MPI_Wtime() - start_time;
    double global_time;
    MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    printf("Rank %d: Time: %lf\n", rank, local_time);
    if(rank == 0)
    {
        printf("Time: %lf\n", global_time);
        print_i_vec(global_bins, b);
    }
    MPI_Finalize();
    return 0;
}