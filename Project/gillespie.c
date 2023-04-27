#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include<unistd.h>
#include "prop.h"

#define m 15
#define n 7
#define T 100
#define R 15

int main(int argc, char *argv[]){

    if(argc != 2){
        printf("Usage %s N\n", argv[0]);
        return -1;
    }
    const int N = atoi(argv[1]);
    int rank, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);    
    
    if(N%num_proc){
        if(rank == 0){
            printf("ERROR N%%p != 0\n");
            return -2;
        }
    }
    const int local_N = N/num_proc;
    
    int X[n*local_N];
    const int P[m*n] ={1, 0, 0, 0, 0, 0, 0,
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
    double w[15];
    srand(time(NULL));
    for(int epoch = 0; epoch < local_N; epoch++){
        int x[n] = {900, 900, 30, 330, 50, 270, 20};
        double t = 0;
        while(t<T){
            prop(x, w);
            double r1 = (double)rand()/RAND_MAX;
            double r2 = (double)rand()/RAND_MAX;
            double tau = -log(r1)/w[14];
            t += tau;
            if(t>T){
                break;
            }
            int j = 0;
            double sum = 0;
            while(sum<r2*w[14]){
                sum += w[j];
                j++;
            }
            j--;
            for(int i = 0; i<n; i++){
                x[i] += P[j*n+i];
            }
        }
        for(int i = 0; i<n; i++){
            X[epoch*n+i] = x[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank == 0){
        int X_all[N*n];
        MPI_Gather(X, n*local_N, MPI_INT, X_all, n*local_N, MPI_INT, 0, MPI_COMM_WORLD);
        FILE *fp = fopen("gillespie.dat", "w");
        for(int i = 0; i<N; i++){
            for(int j = 0; j<n; j++){
                fprintf(fp, "%d ", X_all[i*n+j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }else{
        MPI_Gather(X, n*local_N, MPI_INT, NULL, n*local_N, MPI_INT, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
