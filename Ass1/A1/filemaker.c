#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define func(X) sin(X)

int main(int argc, char *argv[]){

    if(argc != 5){
        printf("Usage %s N a b file_ptr\n", argv[0]);
        printf("Where:\nN is the number of datapoints\n");
        printf("a, b is the function range\n");
        printf("file_ptr is a pointer to file to write to\n");
    }

    const double N = atof(argv[1]);
    const double a = atof(argv[2]);
    const double b = atof(argv[3]);
    char *output_file = argv[4];

    const double dx = (b-a)/(N-1);
    printf("dx: %lf\n", dx);

    double *x = (double *)malloc(N * sizeof(double));
    printf("x: ");
    for(int i = 0; i < N; i++){
        x[i] = i*dx;
        printf("%lf, ", x[i]);
    }
    printf("\n");
    FILE *file;
	if (NULL == (file = fopen(output_file, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
    fprintf(file, "%.4f ", N);

    for (int i = 0; i < N; i++) {
        double value = func(x[i]);
		if (0 > fprintf(file, "%.4f ", value)) {
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