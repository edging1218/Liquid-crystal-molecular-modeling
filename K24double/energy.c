#include "finite.h"

void free_energy(){
	double p_en_ldg = 0;
	double p_en_surf[2] = {0};
        double p_en_el[5] = {0};

	p_en_ldg = energy_ldg();                                                                      
	energy_el(p_en_el);                                                                           
	energy_surf(p_en_surf);                                                                    
	MPI_Barrier(MPI_COMM_WORLD);                                                                  
	MPI_Win_fence(0, win);                                                                        

	//Sum up the energy fraction in different processors and calculate dE
	MPI_Allreduce(&p_en_ldg, &en_ldg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                 
	MPI_Allreduce(p_en_el, en_el, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                     
	MPI_Allreduce(&p_en_surf, &en_surf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);               
//	printf("%d:     %lf\t%lf\t%lf\t%lf.\n", myid, en_ldg, en_el[0], en_el[1], en_surf);   
	en_tot = en_ldg + en_el[0] + en_el[1] + en_el[2] + en_el[3] + en_el[4] + en_surf[0] + en_surf[1];                                  
	dE = en_tot - el_old;                                                             
        el_old = en_tot;
}	

//Landau de-Gennes energy for bulk points.
double energy_ldg(){
	int i, n;
	double ans = 0;
	double traceqq = 0;
	double Qin[6] = {0};
	for (i = 0; i < length; i ++){
		for (n = 0; n < 6; n ++)	Qin[n] = q[i * 6 + n];
		if(sign[i] == 0 || sign[i] == 1){
			traceqq = trqq(Qin);
			ans += 0.5 * (1 - U / 3) * traceqq - U / 3 * trqqq(Qin) + U * 0.25 * traceqq * traceqq;
		}
	}
	return dV * ans;
}

//Elastic energy for bulk points
void energy_el(double* ans){
	int i, n, j, k, l;
	double dQ[3][6];
	double Qin[6] = {0};
	double vec[3] = {0};
	int ref = length * myid;
	double third = 1.0 / 3;
	int xm, xp, ym, yp, zm, zp;
	for (i = 0; i < 3; i ++){
		for(j = 0; j < 6; j ++){
			dQ[i][j] = 0;
		}
	}
	for(i = 0; i < 5; i ++)	ans[i] = 0;

	for (i = 0; i < length; i ++){
		if(sign[i] == 0 || sign[i] == 1){
			for(n = 0; n < 6; n ++)	Qin[n] = q[i * 6 + n];
			xm = neigb[i * 6 + 0] - ref;
			xp = neigb[i * 6 + 1] - ref;
			ym = neigb[i * 6 + 2] - ref;
			yp = neigb[i * 6 + 3] - ref;
			zm = neigb[i * 6 + 4] - ref;
			zp = neigb[i * 6 + 5] - ref;
			for (n = 0; n < 6; n ++) {
				//dQ is the first derivative with second order approximation
				//first index for direction: 0-x; 1-y; 2-z;
				//second index for qtensor index;
				dQ[0][n] = (q[xp * 6 + n] - q[xm * 6 + n]) * 0.5 * idx;
				dQ[1][n] = (q[yp * 6 + n] - q[ym * 6 + n]) * 0.5 * idy;
				dQ[2][n] = (q[zp * 6 + n] - q[zm * 6 + n]) * 0.5 * idz;
			}
			ans[0] += trqq(dQ[0])+trqq(dQ[1])+trqq(dQ[2]);
			if (L2 != 0){
				vec[0] = dQ[0][0];
				vec[1] = dQ[1][1];
				vec[2] = dQ[2][2];
				ans[1] += matr_mult(vec);
				vec[0] = dQ[0][1];
				vec[1] = dQ[1][3];
				vec[2] = dQ[2][4];
				ans[1] += matr_mult(vec);
				vec[0] = dQ[0][2];
				vec[1] = dQ[1][4];
				vec[2] = dQ[2][5];
				ans[1] += matr_mult(vec);
			}
			if (L3 != 0){
				ans[2] += Qin[0] * trqq(dQ[0]) + Qin[3] * trqq(dQ[1]) + Qin[5] * trqq(dQ[2]) \
				        + 2 * Qin[1] * q_mult(dQ[0], dQ[1]) \
					+ 2 * Qin[2] * q_mult(dQ[0], dQ[2]) \ 
					+ 2 * Qin[4] * q_mult(dQ[1], dQ[2]); 
			}
			if (L4 != 0){
				ans[3] += dQ[0][0] * dQ[0][0] + dQ[1][1] * dQ[1][1] + dQ[2][2] * dQ[2][2] \
					+ dQ[0][1] * dQ[0][1] + dQ[1][3] * dQ[1][3] + dQ[2][4] * dQ[2][4] \
					+ dQ[0][2] * dQ[0][2] + dQ[1][4] * dQ[1][4] + dQ[2][5] * dQ[2][5] \
					+ 2 * (dQ[0][1] * dQ[1][0] + dQ[0][2] * dQ[2][0] + dQ[1][2] * dQ[2][1] \
					+ dQ[0][3] * dQ[1][1] + dQ[0][4] * dQ[2][1] + dQ[1][4] * dQ[2][3] \
					+ dQ[1][2] * dQ[0][4] + dQ[0][5] * dQ[2][2] + dQ[1][5] * dQ[2][4]);
			}
			if(chiral == 1){
				//Chiral elastic energy
				ans[4] +=     Qin[0] * dQ[1][2] + Qin[1] * dQ[1][4] + Qin[2] * dQ[1][5] \
					    + Qin[1] * dQ[2][0] + Qin[3] * dQ[2][1] + Qin[4] * dQ[2][2] \
					    + Qin[2] * dQ[0][1] + Qin[4] * dQ[0][3] + Qin[5] * dQ[0][4] \
					    - Qin[0] * dQ[2][1] - Qin[1] * dQ[2][3] - Qin[2] * dQ[2][4] \
					    - Qin[2] * dQ[1][0] - Qin[4] * dQ[1][1] - Qin[5] * dQ[1][2] \
					    - Qin[1] * dQ[0][2] - Qin[3] * dQ[0][4] - Qin[4] * dQ[0][5];
			}
		}
	}
	ans[0] *= 0.5 * dV * L1;
	ans[1] *= 0.5 * dV * L2;
	ans[2] *= 0.5 * dV * L3;
	ans[3] *= 0.5 * dV * L4;
	ans[4] *= dV * chiral * 2 * (L1 - S * third * L3) * qch;
}

void energy_surf(double* ans){
	int i, n, nb;
	double Qdiff[6] = {0};
	double Qin[6] = {0};
	double loc_nu[3] = {0};
	int degen = 0, inf = 1;
	double Wstr = 0;
	double dA = 0;
	bool npboundary = true;
	nb = 0;
	for (i = 0; i < length; i++){
		if(sign[i] >= 2 && sign[i] <= 8){
			//for channel boundary
			if(sign[i] == 2 || sign[i] == 3){
				degen = degenerate;
				inf = infinite;
				Wstr = W;
				dA = dAdrop;
				npboundary = false;
			}
			//for nanoparticle boundary
			else if(sign[i] == 4 || sign[i] == 5){
				degen = 0;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if(sign[i] == 6 || sign[i] == 7){
				degen = 1;
				inf = 0;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else if(sign[i] == 8){
				degen = 0;
				inf = 1;
				Wstr = Wp;
				dA = dApart;
				npboundary = true;
			}
			else{
				printf("Error in energy_surf.\n");
			}
			if(Wstr != 0 && inf != 1){
				//printf("Test. sign[i] = %d\n", sign[i]);
				if(degen == 1){
					for(n = 0; n < 6; n ++)	Qin[n] = q[i * 6 + n];
					for(n = 0; n < 3; n ++)	loc_nu[n] = nu_p[nb * 3 + n];
					en_degen(Qin, loc_nu, Qdiff);
					if(npboundary){
						ans[1] += Wstr * trqq(Qdiff) * dApart;
					}	
					else{
						ans[0] += Wstr * trqq(Qdiff) * dAdrop;
					}
				}
				else if(degen == 0 && inf == 0){
					for(n = 0; n < 6; n ++){
						Qdiff[n] = q[i * 6 + n] - qo_p[nb * 6 + n];
					}
					if(npboundary){
						ans[1] += Wstr * trqq(Qdiff) * dApart;
					}	
					else{
						ans[0] += Wstr * trqq(Qdiff) * dAdrop;
					}
				}
			}
			nb ++;
		}
	}
}

void en_degen(double* Qin, double* loc_nu, double* Qdiff){
	int i, n, j, l, m;
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	double third = 1.0 / 3;
	Qtemp[0][0] = Qin[0] + third * S;
	Qtemp[0][1] = Qtemp[1][0] = Qin[1];
	Qtemp[0][2] = Qtemp[2][0] = Qin[2];
	Qtemp[1][1] = Qin[3] + third * S;
	Qtemp[1][2] = Qtemp[2][1] = Qin[4];
	Qtemp[2][2] = Qin[5] + third * S;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			if(i == j) ptemp[i][j] = 1 - loc_nu[i] * loc_nu[j];
			else ptemp[i][j] = - loc_nu[i] * loc_nu[j];
		}
	}
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			Qp[i][j] = 0;
			for(l = 0; l < 3; l++){
				for(m = 0; m < 3; m++){
					Qp[i][j] += ptemp[i][l]*Qtemp[l][m]*ptemp[m][j];
				}
			}
		}
	}
	Qdiff[0] = Qtemp[0][0] - Qp[0][0];
	Qdiff[1] = Qtemp[0][1] - Qp[0][1];
	Qdiff[2] = Qtemp[0][2] - Qp[0][2];
	Qdiff[3] = Qtemp[1][1] - Qp[1][1];
	Qdiff[4] = Qtemp[1][2] - Qp[1][2];
	Qdiff[5] = Qtemp[2][2] - Qp[2][2];
}
