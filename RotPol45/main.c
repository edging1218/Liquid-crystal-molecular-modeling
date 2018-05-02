#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include "mkl.h"
bool norm_v(float* vec){
	float mod = 0;
	mod = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]); 
	if (mod == 0){
		printf("Zero vector!\n");
		return false;
	}
	else{
		vec[0] = vec[0] / mod;
		vec[1] = vec[1] / mod;
		vec[2] = vec[2] / mod;
	}
	return true;
}
bool checkbound(float* vec, int N){
	int i = 0;
	for(i = 0; i < 3; i ++){
		if(vec[i] < -0.5 || vec[i] > N + 0.5){
			return false;
		}
	}
	return true;
}
void printVec(float* vec){
	int i = 0;
	for(i = 0; i < 3; i ++){
		printf("%f\t", vec[i]);
	}
	printf("\n");
}
void initRot(float A[3][3], float *u, float t){
	t = t * M_PI / 180.0;
	A[0][0] =  cos(t) + u[0] * u[0] * (1 - cos(t));
	A[1][1] =  cos(t) + u[1] * u[1] * (1 - cos(t));
	A[2][2] =  cos(t) + u[2] * u[2] * (1 - cos(t));
	A[0][1] =  u[0] * u[1] * (1 - cos(t)) - u[2] * sin(t);
	A[0][2] =  u[0] * u[2] * (1 - cos(t)) + u[1] * sin(t);
	A[1][0] =  u[1] * u[0] * (1 - cos(t)) + u[2] * sin(t);
	A[1][2] =  u[1] * u[2] * (1 - cos(t)) - u[0] * sin(t);
	A[2][0] =  u[2] * u[0] * (1 - cos(t)) - u[1] * sin(t);
	A[2][1] =  u[2] * u[1] * (1 - cos(t)) + u[0] * sin(t);
	int i, j;
	printf("\nThe rotation matrix is:\n");
	for(i = 0; i < 3; i ++){
		for(j = 0; j < 3; j ++){
			printf("%f\t", A[i][j]);
		}
		printf("\n");
	}
	
}
void Rotate(float A[3][3], float *a, float *b){
	int i, j;
	for(i = 0; i < 3; i ++){
		b[i] = 0;
		for(j = 0; j < 3; j ++){
			b[i] += A[i][j] * a[j];
		}
	}
}
void mult22(double complex *A, double complex *B, double complex *C){
	C[0] = A[0] * B[0] + A[1] * B[2];
	C[1] = A[0] * B[1] + A[1] * B[3];
	C[2] = A[2] * B[0] + A[3] * B[2];
	C[3] = A[2] * B[1] + A[3] * B[3];
	return;
}
void copy22(double complex* A, double complex* B){
	int i;
	for(i = 0; i < 4; i ++){
		A[i] = B[i];
	}
	return;
}
void init22(double complex* A){
	A[0] = 1;
	A[1] = 0;
	A[2] = 0;
	A[3] = 1;
	return;
}
double intensity(double complex* A, double angle){
	angle = angle * M_PI / 180;
	return cabs(- A[0] * sin(angle) * cos(angle) - A[1] * sin(angle) * sin(angle) + A[2] * cos(angle) * cos(angle) + A[3] * sin(angle) * cos(angle));
}

int main(int argc, char *argv[]){
	printf("Valid input should be ./pol ux uy uz theta lambda angle.\n");	
	if(argc != 1 && argc != 2  && argc != 3&& argc != 5 && argc != 6){
		printf("Invalid input.\n");
		return 1;
	}
	float a[9] = {0};
	float v[3] = {0};
	int info;
	int i, n, j, k, l;
	int points;
	int Nx, Ny, Nz;
	int x, y, z;
	double Q[6] = {0};
	double angle = 90;

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
	float* dir;
	dir = (float*)malloc(points * 3 * sizeof(float));

	FILE* grid;
	grid = fopen("grid.bin", "rb");
	if(grid == (FILE*)NULL){
		printf("File grid.bin not found.\n");
		return 1;
	}
	fread(indx, sizeof(int), points, grid);
	fclose(grid);
	//read in the director information
	FILE* qtensor;
	qtensor = fopen("Qtensor.bin", "rb");
	if(qtensor == (FILE*)NULL){
		printf("File Qtensor.bin not found.\n");
		return 1;
	}
	for(i = 0; i < points; i++){
		if(indx[i] == 0 || indx[i] == 1){
			fread(Q, sizeof(double), 6, qtensor);	
			a[0] = (float)Q[0];
			a[3] = (float)Q[1];
			a[6] = (float)Q[2];
			a[4] = (float)Q[3];
			a[7] = (float)Q[4];
			a[8] = (float)Q[5];
			info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, v);
			if (info > 0) {
				printf("Error in eigenvalue routine!\n");
				return 1;
			}
			dir[i * 3] = a[6];
			dir[i * 3 + 1] = a[7];
			dir[i * 3 + 2] = a[8];
		}
		else{
			dir[i * 3] = 0;
			dir[i * 3 + 1] = 0;
			dir[i * 3 + 2] = 0;
		}
	}
	fclose(qtensor);
	
	if(argc == 5 || argc == 6){
		float u[3] = {0};
		u[0] = atof(argv[1]);
		u[1] = atof(argv[2]);
		u[2] = atof(argv[3]);
		if(!norm_v(u)){
			return 1;
		}
		float theta = atof(argv[4]);
		float rotMatrix[3][3] = {0};
		float xyz[3] = {0};
		float rotxyz[3] = {0};
		int lRot, rx, ry, rz;
		float* dirRot;
		dirRot = (float*)malloc(points * 3 * sizeof(float));
		initRot(rotMatrix, u, theta);
		printf("\nConfiguration rotated with u %f, %f, %f and angle %f degree.\n\n", u[0], u[1], u[2], theta);
		l = 0;
		rx = floor(Nx * 0.5);
		ry = floor(Ny * 0.5);
		rz = floor(Nz * 0.5);
		for(k = 0; k < Nz; k++){
			for(j = 0; j < Ny; j++){
				for(i = 0; i < Nx; i++){
					if(indx[l] == 0 || indx[l] == 1){
						xyz[0] = i - rx;
						xyz[1] = j - ry;
						xyz[2] = k - rz;
						Rotate(rotMatrix, xyz, rotxyz);
						rotxyz[0] += rx;
						rotxyz[1] += ry;
						rotxyz[2] += rz;
						if(!checkbound(rotxyz, Nx)){
							printf("Error in rotation.\n");
							printVec(xyz);
							printVec(rotxyz);
						}
						lRot = lrint(round(rotxyz[0])) + lrint(round(rotxyz[1])) * Nx + lrint(round(rotxyz[2])) * Nx * Ny;
						Rotate(rotMatrix, &dir[lRot * 3], &dirRot[l * 3]);
					}
					else{
						dirRot[l * 3] = 0;
						dirRot[l * 3 + 1] = 0;
						dirRot[l * 3 + 2] = 0;
					}
					l ++;
				}
			}
		}
		for(l = 0; l < points * 3; l ++){
			dir[l] = dirRot[l];
		}
		free(dirRot);
		
	}

	FILE* pol;
	double alphap, alpham, alpha;
	double phio, phie, nega;
	double gamma;
	double no, ne, dmesh, lambda;
	no = 1.5;
	ne = 1.7;
	dmesh = 1;
	lambda = 20;
	if(argc == 2 || argc == 3){
		lambda = strtod(argv[1], NULL);
	} 
	else if(argc == 6){
		lambda = strtod(argv[5], NULL);
	}
	if(argc == 3){
		lambda = strtod(argv[1], NULL);
		angle = atof(argv[2]);
	} 
	double complex Pold[4], Pnew[4];
	double complex SR[4];
	pol = fopen("polz.out", "w");
	for(i = 0; i < Nx; i++){
		for(j = 0; j < Ny; j++){
	//		alpham = 0;
			init22(Pold);
			for(k = 0; k < Nz; k++){
				l = i + j * Nx + k * Nx * Ny;
				if(indx[l] == 0){
					gamma = acos(fabs(dir[l * 3 + 2]));	
					alpha = atan2(dir[l * 3 + 1], dir[l * 3]);	
				//	alphap = atan2(dir[l * 3 + 1], dir[l * 3]);	
				//	alpha = alphap - alpham;
				//	alpha = alphap;
					alpham = alphap;
					phio = 2 * M_PI * no * dmesh / lambda;
					nega = no * ne / sqrt(no * no * sin(gamma) * sin(gamma) + ne * ne * cos(gamma) * cos(gamma));
					phie =  2 * M_PI * nega * dmesh / lambda;
/*
					SR[0] = cos(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					SR[1] = -sin(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					SR[2] = sin(alpha) * (cos(phie) + I * sin(phie));
					SR[3] = cos(alpha) * (cos(phie) + I * sin(phie));
*/
					SR[0] = cos(alpha) * cos(alpha) * (cos(phie) + I * sin(phie)) + sin(alpha) * sin(alpha) * (cos(phio) + I * sin(phio));
					SR[1] = sin(alpha) * cos(alpha) * (cos(phie) + I * sin(phie) - cos(phio) - I * sin(phio));
					SR[2] = sin(alpha) * cos(alpha) * (cos(phie) + I * sin(phie) - cos(phio) - I * sin(phio));
					SR[3] = sin(alpha) * sin(alpha) * (cos(phie) + I * sin(phie)) + cos(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					mult22(SR, Pold, Pnew);
					copy22(Pold, Pnew);
				}
			}
		 	fprintf(pol, "%lf\t", intensity(Pold, angle));
		}
		fprintf(pol, "\n");	
	}
	fclose(pol);
	pol = fopen("poly.out", "w");
	for(i = 0; i < Nx; i++){
		for(k = 0; k < Nz; k++){
		//	alpham = 0;
			init22(Pold);
			for(j = 0; j < Ny; j++){
				l = i + j * Nx + k * Nx * Ny;
				if(indx[l] == 0){
					gamma = acos(fabs(dir[l * 3 + 1]));	
				//	alphap = atan2(dir[l * 3], dir[l * 3 + 2]);	
					alpha = atan2(dir[l * 3], dir[l * 3 + 2]);	
				//	alpha = alphap - alpham;
				//	alpham = alphap;
					phio = 2 * M_PI * no * dmesh / lambda;
					nega = no * ne / sqrt(no * no * sin(gamma) * sin(gamma) + ne * ne * cos(gamma) * cos(gamma));
					phie =  2 * M_PI * nega * dmesh / lambda;
/*
					SR[0] = cos(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					SR[1] = -sin(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					SR[2] = sin(alpha) * (cos(phie) + I * sin(phie));
					SR[3] = cos(alpha) * (cos(phie) + I * sin(phie));
*/
					SR[0] = cos(alpha) * cos(alpha) * (cos(phie) + I * sin(phie)) + sin(alpha) * sin(alpha) * (cos(phio) + I * sin(phio));
					SR[1] = sin(alpha) * cos(alpha) * (cos(phie) + I * sin(phie) - cos(phio) - I * sin(phio));
					SR[2] = sin(alpha) * cos(alpha) * (cos(phie) + I * sin(phie) - cos(phio) - I * sin(phio));
					SR[3] = sin(alpha) * sin(alpha) * (cos(phie) + I * sin(phie)) + cos(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					mult22(SR, Pold, Pnew);
					copy22(Pold, Pnew);
				}
			}
		 	fprintf(pol, "%lf\t", intensity(Pold, angle));
		}
		fprintf(pol, "\n");	
	}
	fclose(pol);
	pol = fopen("polx.out", "w");
	for(j = 0; j < Ny; j++){
		for(k = 0; k < Nz; k++){
		//	alpham = 0;
			init22(Pold);
			for(i = 0; i < Nx; i++){
				l = i + j * Nx + k * Nx * Ny;
				if(indx[l] == 0){
					gamma = acos(fabs(dir[l * 3]));	
				//	alphap = atan2(dir[l * 3 + 2], dir[l * 3 + 1]);	
				//	alpha = alphap - alpham;
				//	alpham = alphap;
					alpha = atan2(dir[l * 3 + 2], dir[l * 3 + 1]);	
					phio = 2 * M_PI * no * dmesh / lambda;
					nega = no * ne / sqrt(no * no * sin(gamma) * sin(gamma) + ne * ne * cos(gamma) * cos(gamma));
					phie =  2 * M_PI * nega * dmesh / lambda;
/*
					SR[0] = cos(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					SR[1] = -sin(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					SR[2] = sin(alpha) * (cos(phie) + I * sin(phie));
					SR[3] = cos(alpha) * (cos(phie) + I * sin(phie));
*/
					SR[0] = cos(alpha) * cos(alpha) * (cos(phie) + I * sin(phie)) + sin(alpha) * sin(alpha) * (cos(phio) + I * sin(phio));
					SR[1] = sin(alpha) * cos(alpha) * (cos(phie) + I * sin(phie) - cos(phio) - I * sin(phio));
					SR[2] = sin(alpha) * cos(alpha) * (cos(phie) + I * sin(phie) - cos(phio) - I * sin(phio));
					SR[3] = sin(alpha) * sin(alpha) * (cos(phie) + I * sin(phie)) + cos(alpha) * cos(alpha) * (cos(phio) + I * sin(phio));
					mult22(SR, Pold, Pnew);
					copy22(Pold, Pnew);
				}
			}
		 	fprintf(pol, "%lf\t", intensity(Pold, angle));
		}
		fprintf(pol, "\n");	
	}
			

	free(indx);
	free(dir);
	fclose(pol);
	return 0;
}

