#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "mkl.h"


int main(int argc, char *argv[]){
	if(argc > 2){
		printf("Valid input is ./iso S_threshold.\n");
	}
	float a[9] = {0};
	float v[3] = {0};
	int flag;
	int i, n, j, k, l;
	int info;
	float S;
	float vector[3] = {0};
	int points;
	int Nx, Ny, Nz;
	int x, y, z;
	double Q[6] = {0};

	FILE* param;
	param = fopen("param.in", "r");
	if(param == (FILE*)NULL){
		printf("File param.in not found.\n");
		return 1;
	}
	fscanf(param, "Nx %d\n", &Nx);
	fscanf(param, "Ny %d\n", &Ny);
	fscanf(param, "Nz %d\n", &Nz);
	fclose(param);
	points = Nx * Ny * Nz;

	int* indx;
	indx = (int*)malloc(points * sizeof(int));

	float Sth = 0.60;
	if(argc == 2){
		Sth = atof(argv[1]);
	}

	FILE* grid;
	grid = fopen("grid.bin", "rb");
	if(grid == (FILE*)NULL){
		printf("File grid.bin not found.\n");
		return 1;
	}
	fread(indx, sizeof(int), points, grid);
	fclose(grid);

	FILE* qtensor;
	qtensor = fopen("Qtensor.bin", "rb");
	if(qtensor == (FILE*)NULL){
		printf("File Qtensor.bin not found.\n");
		return 1;
	}
	
	FILE* file;
	file = fopen("Spos.out", "w");
	for (k = 0; k < Nz; k ++) {
		for (j = 0; j < Ny; j++){
			for (i = 0; i < Nx; i++){
				l =  i + Nx * j + Nx * Ny * k;
				if(indx[l] == 0 || indx[l] == 1){
					fread(Q, sizeof(double), 6, qtensor);
					a[0] = Q[0];
					a[3] = Q[1];
					a[6] = Q[2];
					a[4] = Q[3];
					a[7] = Q[4];
					a[8] = Q[5];
					info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, v);
					if (info > 0) {
						printf("Error in eigenvalue routine!\n");
						return 1;
					}
					S = 0.5 * 3 * v[2];
					if(S < Sth){
						fprintf(file, "%d\t%d\t%d\n", i, j, k);
					}
				}
			}
		}
	}
	fclose(qtensor);
	fclose(file);
	free(indx);
	return 0;
}

