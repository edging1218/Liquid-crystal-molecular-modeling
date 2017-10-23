#include "finite.h"

void free_q(){
	//free the space allocated to Qo,Qold
	if(myid == root){
		free(drop);
		free(boundary);
		free(nboundary);
	}
	free(qn);
	if((degenerate == 0 && infinite == 0) || AnchNInf){
		free(qo_p);
	}
	free(nu_p);
//	free(neigb);
	free(sign);
}


bool norm_v(double* vec){
	//Normalize a non-zero vector
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

double dir2ten(double *vec, int n, double Sin){
	//transform dirctor to Qtensor
	double third = 1 / 3.0;
	switch (n) {
		case 0:
			return Sin * (vec[0] * vec[0] - third);
		case 1:
			return Sin * (vec[0] * vec[1]);
		case 2:
			return Sin * (vec[0] * vec[2]);
		case 3:
			return Sin * (vec[1] * vec[1] - third);
		case 4:
			return Sin * (vec[1] * vec[2]);
		case 5:
			return Sin * (vec[2] * vec[2] - third);
		default:
			printf("Error with dir_to_ten!\n");
			break;
			return 0;
	}
	return 0;
}

double trqqq(double* Q){
	//Calculate QijQjkQki
	double ans = 0;
	ans = Q[0] * Q[0] * Q[0] + Q[3] * Q[3] * Q[3] +Q[5] * Q[5] * Q[5]\
		 + 6 * Q[1] * Q[2] * Q[4] + 3 * Q[0] * (Q[1] * Q[1] + Q[2] * Q[2])\
		 + 3 * Q[3] * (Q[1] * Q[1] + Q[4] * Q[4]) + 3 * Q[5] * (Q[4] * Q[4] + Q[2] * Q[2]);
	return ans;
}

double trqq(double* Q){
	//Calculate Q's Frobenius norm//
	double ans = 0;
	ans = Q[0] * Q[0] + Q[3] * Q[3] + Q[5] * Q[5] + 2 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[4] * Q[4]);
	return ans;
}

double q_mult(double* q1, double* q2){
	//Calculate Q1_ijQ2_ij
	double ans = 0;
	ans = q1[0] * q2[0] + q1[3] * q2[3] + q1[5] * q2[5] + 2 * (q1[1] * q2[1] + q1[2] * q2[2] + q1[4] * q2[4]);
	return ans;
}
double matr_mult(double* vec){
	//Calculate norm of vvT
	double ans = 0;
	ans = vec[0] * vec[0] +vec[1] * vec[1] +vec[2] * vec[2] + 2 * vec[0] * vec[1] + 2 * vec[0] * vec[2] + 2 * vec[1] * vec[2]; 
	return ans;
}

bool checktr(double* Q){
	// Check if Qtensor is still traceless, otherwise throw a warning
	double tr = 0;
	double third =  1.0 / 3.0;
	int n;
	tr = (Q[0] + Q[3] + Q[5]) * third;
	if(tr > 1e-5) {
		printf("%f %f %f %f %f %f\n", Q[0], Q[1], Q[2], Q[3], Q[4], Q[5]);
		Q[0] -= tr;
		Q[3] -= tr;
		Q[5] -= tr;
//		printf("Non-tracelss.\n");
		return false;				
	}
	if(trqq(Q) > 1){
//		for(n = 0; n < 6; n ++){
//			Q[n] /= 1.3;
//		}
		printf("Order parameter exceed 1.\n");
		return false;
	}
	return true;
}
