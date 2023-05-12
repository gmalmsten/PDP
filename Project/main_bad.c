#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
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
            double a0 = 0;
            for(int i = 0; i < R; i++){
                a0+=w[i];
            }
            if(a0<0){
                printf("ERROR a0 < 0\nw: [");
                for(int i = 0; i < R; i++){
                    printf("%lf ", w[i]);
                }
                printf("]\nx: [");
                for(int i = 0; i < n; i++){
                    printf("%d ", x[i]);
                }
                printf("]\n");
                return -3;
            }
            double u1 = (double)rand()/RAND_MAX;
            double u2 = (double)rand()/RAND_MAX;
            printf("a0 %lf u1 %lf, u2 %lf\n",a0, u1, u2);

            double tau = -log(u1)/a0;
            double sum = 0;
            double lim = a0*u2;
            int r = 0;
            while(sum < lim){   //-----------------------------------------------------------//
                sum += w[r];
                r++;
                if(r>=R){
                    printf("Error: r exceeds the bounds of the w array\n");
                    return -1;
                }
                if(sum>a0){
                    printf("ERROR: sum > a0");
                    return -1;
                }
            }
            

            for(int i = 0; i < n; i++){
                printf("r: %d x: %d P[i]: %d ",r, x[i], P[r*n + i]);
                x[i] += P[r*n + i];
                printf("x: %d \n", x[i]);
            }
            t+=tau;
            // printf("%lf\n", t);
        }
        for(int i = 0; i < n; i++){
            X[epoch*n + i] = x[i];
        }
        
    }
    
    for(int i = 0; i < local_N; i++){
        printf("Rank %d row %d [", rank, i);
        for(int j = 0; j < n; j++){
            printf("%d ", X[i*n + j]);
        }
        printf("\n");
    }
    MPI_Finalize();
    return 0;
}