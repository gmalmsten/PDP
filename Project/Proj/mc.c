#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "prop.h"
#include <string.h>

// Define macros
#define ROWS 15
#define COLS 7
#define T 100
#define b 20

#define PRODUCE_OUTPUT


void print_i_vec(int *vector, int lim, FILE * restrict fp){
    /*Print a vector consisting of lim integers*/
    fprintf(fp, "[");
    for(int i = 0; i < lim - 1; i++){
        fprintf(fp, "%d, ", vector[i]);
    }
    fprintf(fp, "%d]\n", vector[lim - 1]);
}

int cmp (const void *num1, const void *num2) {
    /*Comparison function for qsort*/
   return ( *(int*)num1 > *(int*)num2 );
}

void print_i_vec_term(int *vector,int lim)
{
    /*Print a vector consisting of lim integers*/
    printf("[");
    for(int i = 0; i < lim - 1; i++){
        printf("%d, ", vector[i]);
    }
    printf("%d]\n", vector[lim - 1]);
}

int main(int argc, char *argv[]){

    if(argc != 3){
        printf("Usage %s N output_file\n", argv[0]);
        return -1;
    }

    // Arguments 
    const int N = atoi(argv[1]);
    const char *output_file = argv[2];
    
    // Initialize MPI
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
    // Uncomment for new seed each run
    // time_t seed = time(NULL);
    // MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    const int seed = 1;
    srand(seed + rank);
    
    // Initialize variables
    int *results = (int *)malloc(local_N*sizeof(int));
    double sub_times[4] = {0}; 
    double all_sub_times[4*num_proc];

    // Initialize window
    MPI_Win win;
    if(num_proc > 1)
    {
        MPI_Win_create(all_sub_times, 4*num_proc*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        MPI_Win_fence(0, win);
    }

    // Start timer
    double start_time = MPI_Wtime();
    for(int epoch = 0; epoch < local_N; epoch++){
        // Initialize new simulation
        double t = 0;
        int x[COLS] = {900, 900, 30, 330, 50, 270, 20};

        char timing = 0;  // Current sub time to store, 0 - [0,25], 1 - (25, 50], 2 - (50, 75], 3 - (75, 100]
        double sub_start_time = MPI_Wtime();
        while(t<T){
            // Accumulate sub times
            if(timing == 0 && t > T/4){
                sub_times[0] += (MPI_Wtime() - sub_start_time);
                sub_start_time = MPI_Wtime();
                timing = 1;
            }
            if(timing == 1 && t > T/2){
                sub_times[1] += (MPI_Wtime() - sub_start_time);
                sub_start_time = MPI_Wtime();
                timing = 2;
            }
            if(timing == 2 && t > 3*T/4){
                sub_times[2] += (MPI_Wtime() - sub_start_time);
                sub_start_time = MPI_Wtime();
                timing = 3;
            }

            // Compute w 
            double w[ROWS];
            prop(x, w);

            // Compute a0
            double a0 = 0;
            for(int i = 0; i < ROWS; i++){
                a0+=w[i];
            }
        

            // Calculate next time step
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
        sub_times[3] += (MPI_Wtime() - sub_start_time);
        results[epoch] = x[0];
    }


    // Rescale sub_times to mean sub_times
    for(int i = 0; i < 4; i++){
        sub_times[i] /= (double)local_N;
    }


    // Put the sub timings in the root process memory
    if(rank == 0)
    {
        memcpy(&all_sub_times[0], sub_times, 4*sizeof(double));
    }
    else
    {
        MPI_Put(sub_times, 4, MPI_DOUBLE, 0, rank * 4, 4, MPI_DOUBLE, win);
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
    MPI_Request min_max_requests[2];
    MPI_Iallreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, &min_max_requests[0]);
    MPI_Iallreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &min_max_requests[1]);

    // Sort the local results while waiting for reduction to complete
    qsort(results, local_N, sizeof(int), cmp);

    // Wait for reduction to complete
    MPI_Waitall(2, min_max_requests, MPI_STATUSES_IGNORE);

    // Close RMA window (Processes are synched following above blocking call)
    if(num_proc>1)
    MPI_Win_free(&win);

    // Calculate bin size and the local counts in each bin
    int bin_size = (global_max - global_min)/b;
    int bins[b] = {0};

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
    free(results);

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
        FILE *fp;
        fp = fopen(output_file, "w");

        fprintf(fp, "Sub-timings:\n");
        fprintf(fp, "Process\t[0, 25]\t\t(25, 50]\t(50, 75]\t(75, 100]\n");
        for(int p = 0; p<num_proc; p++)
        {
            fprintf(fp, "\t%d\t%lf\t%lf\t%lf\t%lf\n", p, all_sub_times[4*p], all_sub_times[4*p+1], all_sub_times[4*p+2], all_sub_times[4*p+3]);
        }
        fprintf(fp, "\nRange of histogram: [%d, %d]\n", global_min, global_max);
        fprintf(fp, "Bins:\n");
        fprintf(fp, "Bin %d [%d %d]\n", 1, global_min, global_min + bin_size);
        for(int i = 1; i < b - 1; i++)
        {
            fprintf(fp, "Bin %d (%d %d]\n", i+1, global_min + bin_size*i, global_min + bin_size*(i + 1));
        }
        fprintf(fp, "Bin %d (%d %d]\n", 20, global_min + bin_size*19, global_max);
        fprintf(fp, "Counts:\n");
        print_i_vec(global_bins, b, fp);
        fclose(fp);
        #endif

        printf("%lf\n", global_time);
    }


    MPI_Finalize();
    return 0;
}