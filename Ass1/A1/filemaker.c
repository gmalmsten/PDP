#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define func(X) sin(X)
#define PI 3.1415926535897932

int main(int argc, char *argv[]){

    if(argc != 3){
        printf("Usage %s N file_ptr\n", argv[0]);
        printf("Where:\nN is the number of datapoints\n");
        printf("file_ptr is a pointer to file to write to\n");
    }

    const int N = atof(argv[1]);
    char *output_file = argv[2];
    
    // if(p<1){
    //     p=1;
    // }

    const double dx = (2*PI)/(N-1);
    printf("dx: %.4lf\n", dx);

    double *x = (double *)malloc(N * sizeof(double));
    for(int i = 0; i < N; i++){
        x[i] = i*dx;
        
    }
    printf("\n");
    FILE *file;
	if (NULL == (file = fopen(output_file, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
    fprintf(file, "%d ", N);

    for (int i = 0; i < N-1; i++) {
        double value = func(x[i]);
		if (0 > fprintf(file, "%lf ", value)) {
			perror("Couldn't write to output file");
		}
	}
    double value = func(x[N-1]);
    if(0 > fprintf(file, "%lf", value)) {
			perror("Couldn't write to output file");
		}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}

    return 0;
}