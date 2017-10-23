#include "finite.h"

void relax_bulk(){
	/* 
	Relax the bulk points with resepct to the phase energy and elastic energy
	*/
	int i, j, n;
	double third = 1.0 / 3.0;

	//for hldg
	double QQ[6] = {0};
	double Qin[6] = {0};
	double Qldg[6] = {0};
	double delta[6] = {1, 0, 0, 1, 0, 1};
	double trQQ = 0;

	//for hel
	double trace = 0;
	double Qelas[6] = {0};
	double Qelas2[6] = {0};
	double Qelas3[6] = {0};
	double ddQ[6][6];
	
	int ref = length * myid;
	
	int xm, xp, ym, yp, zm, zp;
	for (i = 0; i < 6; i ++){
		for(j = 0; j < 6; j ++){
			ddQ[i][j] = 0;
		}
	}
	//for hch
	double Qch[6] ={0};
	double dQ[3][6];
	for (i = 0; i < 3; i ++){
		for(j = 0; j < 6; j ++){
			dQ[i][j] = 0;
		}
	}

	for (i = 0; i < length; i++){
		if(sign[i] == 0 || sign[i] == 1){
			//hldg
			for (n = 0; n < 6; n++) Qin[n] = q[i * 6 + n];
			QQ[0] = Qin[0]*Qin[0]+Qin[1]*Qin[1]+Qin[2]*Qin[2];
			QQ[1] = Qin[0]*Qin[1]+Qin[1]*Qin[3]+Qin[2]*Qin[4];
			QQ[2] = Qin[0]*Qin[2]+Qin[1]*Qin[4]+Qin[2]*Qin[5];
			QQ[3] = Qin[1]*Qin[1]+Qin[3]*Qin[3]+Qin[4]*Qin[4];
			QQ[4] = Qin[1]*Qin[2]+Qin[3]*Qin[4]+Qin[4]*Qin[5];
			QQ[5] = Qin[2]*Qin[2]+Qin[4]*Qin[4]+Qin[5]*Qin[5];
			trQQ = trqq(Qin);
			for (n = 0; n < 6; n++) {
				Qldg[n] = (1-U*third)*Qin[n]-U*(QQ[n]-trQQ*(Qin[n]+delta[n]*third));
			}
		//	if(!checktr2(Qldg)){
		//		printf("Traceless in bulk point. Qldg\n");
		//	}
			xm = neigb[i * 6 + 0] - ref;
			xp = neigb[i * 6 + 1] - ref;
			ym = neigb[i * 6 + 2] - ref;
			yp = neigb[i * 6 + 3] - ref;
			zm = neigb[i * 6 + 4] - ref;
			zp = neigb[i * 6 + 5] - ref;
			for (n = 0; n < 6; n++) {
				//ddQ is de second derivative of qtensor:
				//first index, 0:xx; 1:xy; 2:xz; 3:yy; 4:yz; 5:zz;
				//second index, for qtensor index; 
				ddQ[0][n] = (q[xp * 6 + n]+q[xm * 6 + n]-2*Qin[n])*iddx;
				ddQ[3][n] = (q[yp * 6 + n]+q[ym * 6 + n]-2*Qin[n])*iddy;
				ddQ[5][n] = (q[zp * 6 + n]+q[zm * 6 + n]-2*Qin[n])*iddz;
				Qelas[n] = ddQ[0][n] + ddQ[3][n] + ddQ[5][n];
			}
			if(chiral == 1 || L3 != 0){
				//dQ first derivative: detail see energy.c
				for (n = 0; n < 6; n ++) {
					dQ[0][n] = (q[xp * 6 + n] - q[xm * 6 + n]) * 0.5 * idx;
					dQ[1][n] = (q[yp * 6 + n] - q[ym * 6 + n]) * 0.5 * idy;
					dQ[2][n] = (q[zp * 6 + n] - q[zm * 6 + n]) * 0.5 * idz;
				}
			}
			if((L2 + L4) != 0 || L3 != 0){
				if(sign[i] == 0){
					//neighbor xpyp is the yp neighbor of xp: neighbor[xp * 6 + 3]; same definition for other points
					for (n = 0; n < 6; n ++) {
						ddQ[1][n] = (q[(neigb[xp * 6 + 3] - ref) * 6 + n] + q[(neigb[xm * 6 + 2] - ref) * 6 + n] - q[(neigb[xm * 6 + 3] - ref) * 6 + n] - q[(neigb[xp * 6 + 2] - ref) * 6 + n]) * idx * idy * 0.25;
						ddQ[2][n] = (q[(neigb[xp * 6 + 5] - ref) * 6 + n] + q[(neigb[xm * 6 + 4] - ref) * 6 + n] - q[(neigb[xm * 6 + 5] - ref) * 6 + n] - q[(neigb[xp * 6 + 4] - ref) * 6 + n]) * idx * idz * 0.25;
						ddQ[4][n] = (q[(neigb[yp * 6 + 5] - ref) * 6 + n] + q[(neigb[ym * 6 + 4] - ref) * 6 + n] - q[(neigb[ym * 6 + 5] - ref) * 6 + n] - q[(neigb[yp * 6 + 4] - ref) * 6 + n]) * idy * idz * 0.25;
					}
				}
				else{
					for (n = 0; n < 6; n ++) {
						ddQ[1][n] = 0;
						ddQ[2][n] = 0;
						ddQ[4][n] = 0;
					}
				}
			}
			if((L2 + L4) != 0){
				Qelas2[0] = ddQ[0][0] + ddQ[1][1] + ddQ[2][2];
				Qelas2[1] = 0.5 * (ddQ[1][0] + ddQ[0][1] + ddQ[1][3] + ddQ[3][1] + ddQ[2][4] + ddQ[4][2]);
				Qelas2[2] = 0.5 * (ddQ[2][0] + ddQ[0][2] + ddQ[1][4] + ddQ[4][1] + ddQ[2][5] + ddQ[5][2]);
				Qelas2[3] = ddQ[1][1] + ddQ[3][3] + ddQ[4][4];
				Qelas2[4] = 0.5 * (ddQ[1][2] + ddQ[2][1] + ddQ[5][4] + ddQ[4][5] + ddQ[3][4] + ddQ[4][3]);
				Qelas2[5] = ddQ[2][2] + ddQ[4][4] + ddQ[5][5];
				trace = (Qelas2[0] + Qelas2[3] + Qelas2[5]) * third;
				Qelas2[0] -= trace;
				Qelas2[3] -= trace;
				Qelas2[5] -= trace;
			}
			if(L3 != 0){
				Qelas3[0] = - 0.5 * trqq(dQ[0]);
				Qelas3[1] = - 0.5 * (q_mult(dQ[0], dQ[1]));
				Qelas3[2] = - 0.5 * (q_mult(dQ[0], dQ[2]));
				Qelas3[3] = - 0.5 * trqq(dQ[1]);
				Qelas3[4] = - 0.5 * (q_mult(dQ[1], dQ[2]));
				Qelas3[5] = - 0.5 * trqq(dQ[2]);
				for (n = 0; n < 6; n++){
					Qelas3[n] += Qin[0] * ddQ[0][n] + Qin[3] * ddQ[3][n] + Qin[5] * ddQ[5][n] + 2 * (Qin[1] * ddQ[1][n] + Qin[2] * ddQ[2][n] + Qin[4] * ddQ[4][n]);
					Qelas3[n] += dQ[0][n] * (dQ[0][0] + dQ[1][1] + dQ[2][2]) + dQ[1][n] * (dQ[0][1] + dQ[1][3] + dQ[2][4]) + dQ[2][n] * (dQ[0][2] + dQ[1][4] + dQ[2][5]);
				}		
				trace = (Qelas3[0] + Qelas3[3] + Qelas3[5]) * third;
				Qelas3[0] -= trace;
				Qelas3[3] -= trace;
				Qelas3[5] -= trace;
			}
			if(chiral == 1){
				Qch[0] = 2 * (dQ[1][2] - dQ[2][1]);
				Qch[3] = 2 * (dQ[2][1] - dQ[0][4]);
				Qch[5] = 2 * (dQ[0][4] - dQ[1][2]);
				Qch[1] = dQ[1][4] - dQ[2][3] + dQ[2][0] - dQ[0][2];
				Qch[2] = dQ[1][5] - dQ[2][4] + dQ[0][1] - dQ[1][0];
				Qch[4] = dQ[2][2] - dQ[0][5] + dQ[0][3] - dQ[1][1];
			}
			for (n = 0; n < 6; n++){
				qn[i * 6 + n] = Qin[n] + dt*(- Qldg[n] + L1 * Qelas[n] + (L2 + L4) * Qelas2[n] + L3 * Qelas3[n] - 2 * chiral * qch * (L1 - S * third * L3)* Qch[n]);
			}
		}
	}

	//Wait until all nodes are updated, populat q with qnew
	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Win_fence(0, win);
	for(i = 0; i < length; i ++){	
		if(sign[i] == 0 || sign[i] == 1){
			for (n = 0; n < 6; n++) {
				q[i * 6 + n] = qn[i * 6 + n];
			}
		}
	}
}

void relax_surf(){
	/* 
	Relax the bulk points with resepct to the phase energy, elastic energy and surface energy
	*/
	double dQ[3][6];
	int i, j, n, nb;
	double third = 1.0 / 3.0;
	int degen, inf;
	double Wstr;
	int ref = length * myid;

	//for hel
	double Qelas[6] = {0};
	double Qelas2[6] = {0};
	double Qelas3[6] = {0};
	double Qelas4[6] = {0};
	double trace;
	int xm, xp, ym, yp, zm, zp;
	double temp[3] = {0};

	//for hch
	double Qch[6] ={0};
	for (i = 0; i < 3; i ++){
		for(j = 0; j < 6; j ++){
			dQ[i][j] = 0;
		}
	}

	//for degenerate
	double Qdiff[6] = {0};
	double Qin[6] = {0};
	double loc_nu[3] = {0};

	nb = 0;
	for(i = 0; i < length; i ++){
		if(sign[i] >= 2 && sign[i] <= 8){
			inf = 1;
			//for geo boundary
			if(sign[i] == 2 || sign[i] == 3){
				degen = degenerate;
				inf = infinite;
				Wstr = W;
			}
			//for nanoparticle boundary
			//4 is for NInf
			//6 is for degenerate
			//8 is for Infinite homeotropic
			else if(sign[i] == 4 || sign[i] == 5){
				degen = 0;
				inf = 0;
				Wstr = Wp;
			}
			else if(sign[i] == 6 || sign[i] == 7){
				degen = 1;
				inf = 0;
				Wstr = Wp;
			}
			if(inf == 0){
				for(n = 0; n < 3; n ++)	loc_nu[n] = nu_p[nb * 3 + n];
				for(n = 0; n < 6; n ++)	Qin[n] = q[i * 6 + n];
				xm = neigb[i * 6 + 0] - ref;
				xp = neigb[i * 6 + 1] - ref;
				ym = neigb[i * 6 + 2] - ref;
				yp = neigb[i * 6 + 3] - ref;
				zm = neigb[i * 6 + 4] - ref;
				zp = neigb[i * 6 + 5] - ref;
				if((sign[i] % 2) == 0){
					for (n = 0; n < 6; n++) {
						dQ[0][n] = (-q[xp * 6 + n]+4*q[xm * 6 + n]-3*Qin[n])*0.5*idx;
						dQ[1][n] = (-q[yp * 6 + n]+4*q[ym * 6 + n]-3*Qin[n])*0.5*idy;
						dQ[2][n] = (-q[zp * 6 + n]+4*q[zm * 6 + n]-3*Qin[n])*0.5*idz;
						Qelas[n] = dQ[0][n] * fabs(loc_nu[0]) + dQ[1][n] * fabs(loc_nu[1]) + dQ[2][n] * fabs(loc_nu[2]);
					}
				}
				else if((sign[i] % 2) == 1){
					for (n = 0; n < 6; n++) {
						if((xm + ref) == -1){
							dQ[0][n] = 0;
						}
						else if((xp + ref) == -1){
							dQ[0][n] = (q[xm * 6 + n]-Qin[n])*idx;
						}
						else{
							dQ[0][n] = (-q[xp * 6 + n]+4*q[xm * 6 + n]-3*Qin[n])*0.5*idx;
						}
						if((ym + ref) == -1){
							dQ[1][n] = 0;
						}
						else if((yp + ref) == -1){
							dQ[1][n] = (q[ym * 6 + n]-Qin[n])*idy;
						}
						else{
							dQ[1][n] = (-q[yp * 6 + n]+4*q[ym * 6 + n]-3*Qin[n])*0.5*idy;
						}
						if((zm + ref) == -1){
							dQ[2][n] = 0;
						}
						else if((zp + ref) == -1){
							dQ[2][n] = (q[zm * 6 + n]-Qin[n])*idz;
						}
						else{
							dQ[2][n] = (-q[zp * 6 + n]+4*q[zm * 6 + n]-3*Qin[n])*0.5*idz;
						}
						Qelas[n] = dQ[0][n] * fabs(loc_nu[0]) + dQ[1][n] * fabs(loc_nu[1]) + dQ[2][n] * fabs(loc_nu[2]);
					}
				}
		//		if(!checktr2(Qelas)){
		//			printf("Traceless in surf point. Qelas\n");
		//		}
				if(L2 != 0 || L3 != 0 || L4 != 0){
					for(j = 0; j < 3; j++){
						if(loc_nu[j] < 0){
							for(n = 0; n < 6; n++){
								dQ[j][n] = -dQ[j][n];
							}
						}
					}
				}
				if(L2 != 0){
					temp[0] = dQ[0][0] + dQ[1][1] + dQ[2][2];
					temp[1] = dQ[0][1] + dQ[1][3] + dQ[2][4];
					temp[2] = dQ[0][2] + dQ[1][4] + dQ[2][5];
					trace = (loc_nu[0] * temp[0] + loc_nu[1] * temp[1] + loc_nu[2] * temp[2]) * third;
					Qelas2[0] = loc_nu[0] * temp[0] - trace;
					Qelas2[3] = loc_nu[1] * temp[1] - trace;
					Qelas2[5] = loc_nu[2] * temp[2] - trace;
					Qelas2[1] = 0.5 * (loc_nu[0] * temp[1] + loc_nu[1] * temp[0]);
					Qelas2[2] = 0.5 * (loc_nu[0] * temp[2] + loc_nu[2] * temp[0]);
					Qelas2[4] = 0.5 * (loc_nu[2] * temp[1] + loc_nu[1] * temp[2]);
				}
				if(L3 != 0){
					for(n = 0; n < 6; n++){
						Qelas3[n] = Qin[0] * dQ[0][n] * loc_nu[0] + Qin[3] * dQ[1][n] * loc_nu[1] + Qin[5] * dQ[2][n] * loc_nu[2]\
							    + Qin[1] * (dQ[0][n] * loc_nu[1] + dQ[1][n] * loc_nu[0])\
							    + Qin[2] * (dQ[0][n] * loc_nu[2] + dQ[2][n] * loc_nu[0])\
							    + Qin[4] * (dQ[1][n] * loc_nu[2] + dQ[2][n] * loc_nu[1]);
					}
					trace = (Qelas3[0] + Qelas3[3] + Qelas3[5]) * third;
					Qelas3[0] -= trace;
					Qelas3[3] -= trace;
					Qelas3[5] -= trace;
				}
				if(L4 != 0){
					Qelas4[0] = dQ[0][0] * loc_nu[0] + dQ[0][1] * loc_nu[1] + dQ[0][2] * loc_nu[2];
					Qelas4[3] = dQ[1][1] * loc_nu[0] + dQ[1][3] * loc_nu[1] + dQ[1][4] * loc_nu[2];
					Qelas4[5] = dQ[2][2] * loc_nu[0] + dQ[2][4] * loc_nu[1] + dQ[2][5] * loc_nu[2];
					Qelas4[1] = 0.5 * (loc_nu[0] * (dQ[0][1] + dQ[1][0]) + loc_nu[1] * (dQ[1][1] + dQ[0][3]) + loc_nu[2] * (dQ[1][2] + dQ[0][4]));
					Qelas4[2] = 0.5 * (loc_nu[0] * (dQ[0][2] + dQ[2][0]) + loc_nu[1] * (dQ[2][1] + dQ[0][4]) + loc_nu[2] * (dQ[2][2] + dQ[0][5]));
					Qelas4[4] = 0.5 * (loc_nu[0] * (dQ[2][1] + dQ[1][2]) + loc_nu[1] * (dQ[1][4] + dQ[2][3]) + loc_nu[2] * (dQ[1][5] + dQ[2][4]));
					trace = (Qelas4[0] + Qelas4[3] + Qelas4[5]) * third;
					Qelas4[0] -= trace;
					Qelas4[3] -= trace;
					Qelas4[5] -= trace;
				}
				if(chiral == 1){
					Qch[0] = loc_nu[2] * Qin[1] - loc_nu[1] * Qin[2];
					Qch[3] = loc_nu[0] * Qin[4] - loc_nu[2] * Qin[1];
					Qch[5] = loc_nu[1] * Qin[2] - loc_nu[0] * Qin[4];
					Qch[1] = 0.5 * (loc_nu[2] * Qin[3] - loc_nu[1] * Qin[4] + loc_nu[0] * Qin[2] - loc_nu[2] * Qin[0]);
					Qch[2] = 0.5 * (loc_nu[2] * Qin[4] - loc_nu[1] * Qin[5] + loc_nu[1] * Qin[0] - loc_nu[0] * Qin[1]);
					Qch[4] = 0.5 * (loc_nu[0] * Qin[5] - loc_nu[2] * Qin[2] + loc_nu[1] * Qin[1] - loc_nu[0] * Qin[3]);
				}
			}
			if(degen == 1){
				relax_degen(Qin, loc_nu, Qdiff);
				for(n = 0; n < 6; n++){
					qn[i * 6 + n] = Qin[n] + dt*(L1 * Qelas[n] + L2 * Qelas2[n] + L3 * Qelas3[n] + L4 * Qelas4[n] + chiral * 2 * qch * (L1 - S * third * L3)* Qch[n] - 2 * Wstr * Qdiff[n]);
				}
			}
			else if(degen == 0 && inf == 0){
				for(n = 0; n < 6; n++){
					qn[i * 6 + n] = Qin[n] + dt*(L1 * Qelas[n] + L2 * Qelas2[n] + L3 * Qelas3[n] + L4 * Qelas4[n] + chiral * 2 * qch * (L1 - S * third * L3)* Qch[n] - Wstr* (Qin[n]-qo_p[nb * 6 + n]));
				}
			}
			nb ++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Win_fence(0, win);
	for(i = 0; i < length; i ++){	
		if((sign[i] >= 4 && sign[i] < 8) || ((sign[i] == 2 || sign[i] == 3) && infinite == 0)){
			for (n = 0; n < 6; n++) {
				q[i * 6 + n] = qn[i * 6 + n];
			}
		}
	}
}

void relax_degen(double* Qin, double* loc_nu, double* Qdiff){
	int i, n, j, l, m;
	double Qtemp[3][3];
	double ptemp[3][3];
	double Qp[3][3];
	double third = 1.0 / 3;
	double nuQnu;
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
	nuQnu = 0;
	for(i = 0; i<3; i++) {
		for(j = 0; j<3; j++){
			nuQnu += loc_nu[i]*Qtemp[i][j]*loc_nu[j];
		}
	}
	nuQnu *= third;
	for(n = 0; n < 6; n ++)	Qdiff[n] = 0;
	Qdiff[0] =  Qtemp[0][0]- Qp[0][0] - nuQnu;
	Qdiff[1] =  Qtemp[0][1]- Qp[0][1];
	Qdiff[2] =  Qtemp[0][2]- Qp[0][2];
	Qdiff[3] =  Qtemp[1][1]- Qp[1][1] - nuQnu;
	Qdiff[4] =  Qtemp[1][2]- Qp[1][2];
	Qdiff[5] =  Qtemp[2][2]- Qp[2][2] - nuQnu;
}
