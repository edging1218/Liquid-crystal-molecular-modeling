#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "mkl.h"



int main(int argc, char *argv[]){
	float a[9] = {0};
	float v[3] = {0};
	int flag;
	int i, n, j, k, l;
	int info;
	int points;
	int Nx, Ny, Nz;
	float x, y, z;

	//Read in parameter
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
	int rx = lrint(Nx * 0.5);
	int ry = lrint(Ny * 0.5);
	int rz = lrint(Nz * 0.5);
	
	double Qloc[6] = {0};
	int indx = 0;

	//define output and initialize as zero
	float* S_array_x, *S_array_y, *S_array_z;
	S_array_x = (float*)malloc(Nx * sizeof(float));
	S_array_y = (float*)malloc(Ny * sizeof(float));
	S_array_z = (float*)malloc(Nz * sizeof(float));
	float* TwAn_x, *TwAn_y, *TwAn_z;
	TwAn_x = (float*)malloc(Nx * sizeof(float));
	TwAn_y = (float*)malloc(Ny * sizeof(float));
	TwAn_z = (float*)malloc(Nz * sizeof(float));
	for(i = 0; i < Nx; i ++){
		S_array_x[i] = 0;
		TwAn_x[i] = 0;
	}
	for(i = 0; i < Ny; i ++){
		S_array_y[i] = 0;
		TwAn_y[i] = 0;
	}
	for(i = 0; i < Nz; i ++){
		S_array_z[i] = 0;
		TwAn_z[i] = 0;
	}

	FILE* grid;
	grid = fopen("grid.bin", "rb");
	if(grid == (FILE*)NULL){
		printf("File grid.bin not found.\n");
		return 1;
	}
	FILE* qtensor;
	qtensor = fopen("Qtensor.bin", "rb");
	if(qtensor == (FILE*)NULL){
		printf("File Qtensor.bin not found.\n");
		return 1;
	}

	for (k = 0; k < Nz; k ++) {
		for (j = 0; j < Ny; j++){
			for (i = 0; i < Nx; i++){
				fread(&indx, sizeof(int), 1, grid);
				if(indx == 0 || indx == 1){
					fread(Qloc, sizeof(double), 6, qtensor);
					if(k == rz && j == ry){
						a[0] = (float)Qloc[0];
						a[3] = (float)Qloc[1];
						a[6] = (float)Qloc[2];
						a[4] = (float)Qloc[3];
						a[7] = (float)Qloc[4];
						a[8] = (float)Qloc[5];
						info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, v);
						if (info > 0) {
							printf("Error in eigenvalue routine!\n");
							return 1;
						}
						S_array_x[i] = 0.5 * 3 * v[2];
						TwAn_x[i] = a[8];
					}
					if(i == rx && j == ry){
						a[0] = (float)Qloc[0];
						a[3] = (float)Qloc[1];
						a[6] = (float)Qloc[2];
						a[4] = (float)Qloc[3];
						a[7] = (float)Qloc[4];
						a[8] = (float)Qloc[5];
						info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, v);
						if (info > 0) {
							printf("Error in eigenvalue routine!\n");
							return 1;
						}
						S_array_z[i] = 0.5 * 3 * v[2];
						TwAn_z[i] = a[8];
					}
					if(k == rz && i == rx){
						a[0] = (float)Qloc[0];
						a[3] = (float)Qloc[1];
						a[6] = (float)Qloc[2];
						a[4] = (float)Qloc[3];
						a[7] = (float)Qloc[4];
						a[8] = (float)Qloc[5];
						info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, v);
						if (info > 0) {
							printf("Error in eigenvalue routine!\n");
							return 1;
						}
						S_array_y[i] = 0.5 * 3 * v[2];
						TwAn_y[i] = a[8];
					}
				}

			}
		}
	}
	fclose(qtensor);
	fclose(grid);

	//Write to file
	FILE* file;
	file=fopen("Sx.out","w");
	for(i= 0; i< Nx; i ++){
		fprintf(file, "%d %f %f\n", i, S_array_x[i], fabs(TwAn_x[i]));
	}	
	fclose(file);
	file=fopen("Sy.out","w");
	for(i= 0; i< Ny; i ++){
		fprintf(file, "%d %f %f\n", i, S_array_y[i], fabs(TwAn_y[i]));
	}	
	fclose(file);
	file=fopen("Sz.out","w");
	for(i= 0; i< Nz; i ++){
		fprintf(file, "%d %f %f\n", i, S_array_z[i], fabs(TwAn_z[i]));
	}	
	fclose(file);
	//free 
	free(S_array_x);
	free(S_array_y);
	free(S_array_z);
	free(TwAn_x);
	free(TwAn_y);
	free(TwAn_z);
	return 0;
}

