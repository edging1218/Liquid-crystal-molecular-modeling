#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "mkl.h"
#define PI 3.1415926

bool period[3] = {false, false, false};
int Nx, Ny, Nz;

double dot(double* n1, double* n2){
	return n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
}

int peri(int idx, int dir){
	if(!period[dir]) return idx;
	if(dir == 0){
		if(idx < 0)	
			return idx + Nx;
		else if(idx >= Nx)
			return idx - Nx;
		else
			return idx;
	}
	else if(dir == 1){
		if(idx < 0)	
			return idx + Ny;
		else if(idx >= Ny)
			return idx - Ny;
		else 
			return idx;
	}
	else if(dir == 2){
		if(idx < 0)	
			return idx + Nz;
		else if(idx >= Nz)
			return idx - Nz;
		else 
			return idx;
	}
	else{
		printf("Wrong input for function periodic.\n");
	}
	return 0;
}

void reverse(double* n1){
	n1[0] = -n1[0];
	n1[1] = -n1[1];
	n1[2] = -n1[2];
}

bool norm_v(double* vec){
	double mod = 0;
	mod = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]); 
	if (mod == 0){
		printf("Zero vector!\n");
		return false;
	}
	else{
		vec[0] = vec[0] / mod;
		vec[1] = vec[1] / mod;
		vec[2] = vec[2] / mod;
		return true;
	}
}

int main(int argc, char *argv[]){
	float a[9] = {0};
	float v[3] = {0};
	double Q[6] = {0};
	bool flag = true;
	int i, j, k, l, n;
	int info;
	int points;
	int x, y, z;
	double splay, bend, twist;
	splay = 0;
	bend = 0;
	twist = 0;
	double twist_temp = 0;
	double bend_temp[3] = {0};

	int xm, xp, ym, yp, zm, zp;
	double nxm[3], nxp[3], nym[3], nyp[3], nzm[3], nzp[3], ncur[3];
	double dn[3][3];
	double div;
	double curl[3];
	int bulk = 0;
	//	printf("The threshold for order parameter.\n");
	//	scanf("%lf%*c", &Sth);


	double Lx, Ly, Lz;;
	double W, U, L1, L2;
	int chiral;
	double N;
	FILE* param;
	param = fopen("param.in", "r");
	if(param == (FILE*)NULL){
		printf("File param.in not found.\n");
		return 1;
	}
	fscanf(param, "Nx %d\n", &Nx);
	fscanf(param, "Ny %d\n", &Ny);
	fscanf(param, "Nz %d\n", &Nz);
	fscanf(param, "Lx %lf\n", &Lx);
	fscanf(param, "Ly %lf\n", &Ly);
	fscanf(param, "Lz %lf\n", &Lz);
	fscanf(param, "W %lf\n", &W);
	fscanf(param, "U %lf\n", &U);
	fclose(param);
	points = Nx * Ny * Nz;
	double qch = 0;
//	double qch = PI * N / 2.0 / D;
	double S = 0.25 * (1 + 3 * sqrt(1 - 8.0 / (3 * U)));
	printf("Seq is %lf.\n", S);
//	printf("qch is %lf.\n", qch);
	//	double Sth = S;

	int* indx;
	indx = (int*)malloc(points * sizeof(int));
	double* dir;
	dir = (double*)malloc(points * 3 * sizeof(double));
	double* S_array;
	S_array = (double*)malloc(points * sizeof(double));

//	FILE* file;
//	file = fopen("splay.vtk", "w");
//	if(file == (FILE*)NULL)
//		return 1;
//	fprintf(file, "# vtk DataFile Version 3.0\n");
//	fprintf(file, "vtk output\n");
//	fprintf(file, "ASCII\n");
//	fprintf(file, "DATASET STRUCTURED_GRID\n");
//	fprintf(file, "DIMENSIONS\t %d\t %d\t %d\t\n", Nx, Ny, Nz);
//	fprintf(file, "POINTS\t %d\t float\n", points);

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
	for(i = 0; i < points; i++){
		if(indx[i] == 0 || indx[i] == 1) {
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
			dir[i * 3 + 0] = (double)a[6];
			dir[i * 3 + 1] = (double)a[7];
			dir[i * 3 + 2] = (double)a[8];
			if(!norm_v(&dir[i * 3])){
				printf("Error in director.\n");
				flag = false;	
				break;
			}
			S_array[i] = 0.5 * 3 * v[2];
			bulk ++;
		}
		else{
			dir[i * 3] = 0;
			dir[i * 3 + 1] = 0;
			dir[i * 3 + 2] = 0;
			S_array[i] = 0;
		}
	}
	fclose(qtensor);

	printf("There are %d total points.\n", points);
	printf("There are %d droplet points.\n", bulk);
	//	bulk = 0;
//	fprintf(file, "\n");
//	fprintf(file, "POINT_DATA\t%d\n", points);
//	fprintf(file, "\n");
//	fprintf(file, "SCALARS splay float 1\n");
//	fprintf(file, "LOOKUP_TABLE default\n");
	if(flag){
		for (k = 0; k < Nz; k ++) {
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					l = i + j * Nx + k * Nx * Ny;
					div = 0;
					if(indx[l] == 0 && S_array[l] >= 0.6){
						xm = peri(i - 1, 0) + j * Nx + k * Nx * Ny;
						xp = peri(i + 1, 0) + j * Nx + k * Nx * Ny;
						ym = i + peri(j - 1, 1) * Nx + k * Nx * Ny;
						yp = i + peri(j + 1, 1) * Nx + k * Nx * Ny;
						zm = i + j * Nx + peri(k - 1, 2) * Nx * Ny;
						zp = i + j * Nx + peri(k + 1, 2) * Nx * Ny;
						if(xm >= points || xp >= points || ym >= points || yp >= points || zm >= points || zp >= points){
							printf("%d %d %d\n", i, j, k);
							flag = false;
							break;
						}
						for (n = 0; n < 3; n++) {
							ncur[n] = dir[l * 3 + n];
							nxm[n] = dir[xm * 3 + n];
							nxp[n] = dir[xp * 3 + n];
							nym[n] = dir[ym * 3 + n];
							nyp[n] = dir[yp * 3 + n];
							nzm[n] = dir[zm * 3 + n];
							nzp[n] = dir[zp * 3 + n];
						}	
						if(dot(ncur, nxm) < 0){
							reverse(nxm);
						}
						if(dot(ncur, nxp) < 0){
							reverse(nxp);
						}
						if(dot(ncur, nym) < 0){
							reverse(nym);
						}
						if(dot(ncur, nyp) < 0){
							reverse(nyp);
						}
						if(dot(ncur, nzm) < 0){
							reverse(nzm);
						}
						if(dot(ncur, nzp) < 0){
							reverse(nzp);
						}
						for (n = 0; n < 3; n++) {
							dn[0][n] = (nxp[n] - nxm[n]) * 0.5;
							dn[1][n] = (nyp[n] - nym[n]) * 0.5;
							dn[2][n] = (nzp[n] - nzm[n]) * 0.5;
						}
						div = dn[0][0] + dn[1][1] + dn[2][2];
						splay += div * div; 
						curl[0] = dn[1][2] - dn[2][1];
						curl[1] = dn[2][0] - dn[0][2];
						curl[2] = dn[0][1] - dn[1][0];
						twist_temp = dot(ncur, curl);
						twist += (twist_temp + qch) * (twist_temp + qch);
						bend_temp[0] = ncur[1] * curl[2] - ncur[2] * curl[1];
						bend_temp[1] = ncur[2] * curl[0] - ncur[0] * curl[2];
						bend_temp[2] = ncur[0] * curl[1] - ncur[1] * curl[0];
						bend += dot(bend_temp, bend_temp);
					}
				 //	fprintf(file, "%lf\n", div);
				}	
				if(!flag)	break;
			}
			if(!flag)	break;
		}
	}
//	printf("There are %d bulk points.\n", bulk);
	FILE* elas;
	elas = fopen("elas.out", "w");
	double tot = splay + twist + bend;
	printf("\nN\tq\tSplay\ttwist\tbend\ttot\tr_splay\tr_twist\tr_bend\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", N, qch, splay, twist, bend, tot, splay/tot, twist/tot, bend/tot);
	fprintf(elas, "\nN\tq\tSplay\ttwist\tbend\ttot\tr_splay\tr_twist\tr_bend\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", N, qch, splay, twist, bend, tot, splay/tot, twist/tot, bend/tot);
	fclose(elas);
//	fclose(file);
	free(dir);
	free(indx);
	free(S_array);
	return 0;
}

