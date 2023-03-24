#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define func(X) sin(X)
#define PI 3.1415926535897932

int main(int argc, char *argv[]){

    if(argc != 5){
        printf("Usage %s N n file_ptr\n", argv[0]);
        printf("Where:\nN is the number of datapoints\n");
        printf("n is the number of periods\n");
        printf("file_ptr is a pointer to file to write to\n");
    }

    const double N = atof(argv[1]);
    const double n = atof(argv[2]);
    char *output_file = argv[3];

    const double dx = (2*PI)/(N-1);
    printf("dx: %lf\n", dx);

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
    fprintf(file, "%d ", (int)N);

    for (int i = 0; i < N; i++) {
        double value = func(x[i]);
		if (0 > fprintf(file, "%f ", value)) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}

    return 0;
}