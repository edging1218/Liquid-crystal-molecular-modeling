#include "finite.h"

bool initial_cylinder(){
	int l;
	int nb, nd, np;
	int i, j, k, n, m;
	double dis, x, y, z;
	double x1, y1, z1;
	int rx = lrint(Nx / 2);
	int ry = lrint(Ny / 2);
	double R = rx - 2;
	int xm, xp, ym, yp, zm, zp;
	int nsurf;
	bool *ndrop;
	int *indx;
	double **pos;
	int count1;
	double dir[3] = {0, 0, 1};
	
	double dx = Lx/(Nx-1);
	double dy = Ly/(Ny-1);
	double dz = Lz/(Nz-1);

	idx = 1 / dx;
	idy = 1 / dy;
	idz = 1 / dz;
	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

	bulk = 0;
	surf = 0;
	nsurf = 0;
	tot = Nx * Ny * Nz;
        qch = 0.5 * N / R * M_PI;

	//allocate drop and boundary
	ndrop = (bool*)malloc(tot * sizeof(bool));
	nboundary = (bool*)malloc(tot * sizeof(bool));
	drop = (bool*)malloc(tot * sizeof(bool));
	boundary = (bool*)malloc(tot * sizeof(bool));
	indx = (int*)malloc(tot * sizeof(int));
	for(l = 0; l < tot; l ++){
		drop[l] = false;
		boundary[l] = false;
		ndrop[l] = false;
		nboundary[l] = false;
		indx[l] = -1;
	}

	l = 0;
	//define the droplet 
	for(k = 0; k < Nz; k++){
		for (j = 0; j < Ny; j++){
			for (i = 0; i < Nx; i++){
				x = (i-rx)*dx;
				y = (j-ry)*dy;
				dis = sqrt(x*x+y*y);
				if (dis <= (R + 0.5)){
					drop[l] = true;
					bulk ++;
				}
				l ++;
			}
		}
	}
	droplet = bulk;

	//define boundary
	l = 0;
	for(k = 0; k < Nz; k++){
		for (j = 0; j < Ny; j++){
			for (i = 0; i < Nx; i++){
				if(drop[l]){
					xm = i - 1 + j * Nx + k * Nx * Ny;
					xp = i + 1 + j * Nx + k * Nx * Ny;
					ym = i + (j - 1) * Nx + k * Nx * Ny;
					yp = i + (j + 1) * Nx + k * Nx * Ny;
					if(!drop[xm] || !drop[xp] || !drop[ym] || !drop[yp]){
						boundary[l] = true;
						surf ++;
					}
				}
				l ++;
			}
		}
	}
	bulk -= surf;
	for(l = 0; l < tot; l ++){
		if(boundary[l])		drop[l] = false;
	}


	//Initialize particle if Np > 0
	if(Np != 0){
		pos = (double**)malloc(Np * sizeof(double*));
		for(i = 0; i < Np; i ++){
			pos[i] = (double*)malloc(4 * sizeof(double));
			for(n = 0; n < 4; n ++){
				pos[i][n] = 0;
			}
		}
		if(!read_nppos(pos)) return false;

		l = 0;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(boundary[l] || drop[l]){
						for(m = 0; m < Np; m ++){
							x = i - pos[m][0];
							y = j - pos[m][1];
							z = k - pos[m][2];
							if(z > 0.5 * Nz){
								z -= Nz;
							}
							else if(z < -0.5 * Nz){
								z += Nz;
							}
							dis = x*x+y*y+z*z;
							if (dis <= (Rp - 0.5) * (Rp - 0.5)){
								ndrop[l] = true;
								if(drop[l]){
									bulk --;
									drop[l] = false;
								}
								else{	
									boundary[l] = false;
									surf --;
								}
							}
						}						
					}
					l ++;
				}
			}
		}
		l = 0;
		for(k = 0; k < Nz; k++){
			for (j = 0; j < Ny; j++){
				for (i = 0; i < Nx; i++){
					if(drop[l] || boundary[l]){
						xm = i - 1 + j * Nx + k * Nx * Ny;
						xp = i + 1 + j * Nx + k * Nx * Ny;
						ym = i + (j - 1) * Nx + k * Nx * Ny;
						yp = i + (j + 1) * Nx + k * Nx * Ny;
						zm = i + j * Nx + peri(k - 1, 2) * Nx * Ny;
						zp = i + j * Nx + peri(k + 1, 2) * Nx * Ny;
						if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zm] || ndrop[zp]){
							if(drop[l]){
								nboundary[l] = true;
								drop[l] = false;
								bulk --;
								nsurf ++;	
								surf ++;
							}	
							else if(Wp > W){
								boundary[l] = false;
								nboundary[l] = true;
								nsurf ++;
							}
						}
					}
					l ++;
				}
			}
		}
	}
	//calculate dV and dA for NP on the surf
//	dV = (4.0 / 3 * M_PI * R * R * R - 4.0 / 3 * M_PI * Rp * Rp * Rp * Np * 0.5) / bulk; 
//	dAdrop = (4 * M_PI * R * R - M_PI * Rp * Rp * Np) / (surf -  nsurf);
//	dApart = (2 * M_PI * Rp * Rp * Np) / nsurf ;
//	printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 

	//calculate dV and dA for NP inside droplet
	dV = (M_PI * R * R * Lz - 4.0 / 3 * M_PI * Rp * Rp * Rp * Np) / bulk; 
	dAdrop = (2 * M_PI * R * Lz) / (surf -  nsurf);
	if(Np != 0){
		dApart = (4 * M_PI * Rp * Rp * Np) / nsurf ;
	}
	else{
		dApart = 0;
	}

	printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 
	droplet = bulk + surf;
	printf("\nR is %lf\nDroplet nodes number is %d\nBulk nodes number is %d\nDroplet surface nodes number is %d\nParticle surface nodes number is %d\n", R, droplet, bulk, surf, nsurf); 

	//allocate nu 
	//allocate Qo only for finite homeotropic
	nu = (double*)malloc(surf * 3 * sizeof(double));
	for(i = 0; i < surf * 3; i ++){
		nu[i] = 0;
	}
	AnchNInf = false;
	for(i = 0; i < Np; i++){
		if(pos[i][3] == 1){
			AnchNInf = true;
			break;
		}
	}
	if((degenerate == 0 && infinite == 0) || AnchNInf){
		Qo = (double*)malloc(6 * surf * sizeof(double));
		for(i = 0; i < surf * 6; i ++){
			Qo[i] = 0;
		}
	}

	//allocate qold and neighbor
	//allocate share to define droplet: -1 not defined; 20 bulk; 0-9 droplet boundary; 10 -19 nanoparticle boundary
	length = lrint(droplet / numprocs) + 1;
	share = (int*)malloc(numprocs * length * sizeof(int));
	Qold = (double*)malloc(6 * numprocs * length * sizeof(double));
	neighbor = (int*)malloc(6 * numprocs * length * sizeof(int));
	for(i = 0; i < numprocs * length; i ++){
		share[i] = -1;
	}
	for(i = 0; i < 6 * numprocs * length; i ++){
		Qold[i] = 0;
		neighbor[i] = -1;
	}

	//populate indx array to transformation from 3D to 1D and share array.
	nb = 0;
	nd = 0;
	for(l = 0; l < tot; l++){
		if(drop[l] || boundary[l] || nboundary[l]){
			indx[l] = nd;
			if(drop[l]) share[nd] = 0;
			else if(boundary[l] || nboundary[l]){
				if(boundary[l])	share[nd] = 2;
				else	share[nd] = 4;
				nb ++;
			}
			nd ++;
		}
	}

	if (nd != droplet){
		printf("Problem in initialization of qtensor. nd is %d not equal to droplet %d.\n", nd, droplet);
		return false;
	}

	if (nb != surf){
		printf("Problem in initialization of qtensor. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}

	//Initial configuration of qtensor according to seed 
	if(!conf(pos))	return false;

	//define neighbors and calculate normal vector nu for droplet and particle.
	l = 0;
	nb = 0;
	for(k = 0; k < Nz; k++){
		for (j = 0; j < Ny; j++){
			for (i = 0; i < Nx; i++){
				nd = indx[l];
				if(drop[l]){
					neighbor[nd * 6 + 0] = indx[i - 1 + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 1] = indx[i + 1 + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 2] = indx[i + (j - 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 3] = indx[i + (j + 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 4] = indx[i + j * Nx + peri(k - 1, 2) * Nx * Ny];
					neighbor[nd * 6 + 5] = indx[i + j * Nx + peri(k + 1, 2) * Nx * Ny];
				}
				if(boundary[l] || nboundary[l]){
					if(boundary[l]){
						x = (i-rx)*dx;
						y = (j-ry)*dy;
						dis = sqrt(x*x+y*y);
						//define nu
						if (dis == 0){
							printf("Error in neighbors on boundary.\n");
							return false;
						}
						nu[nb * 3 + 0] = -x;
						nu[nb * 3 + 1] = -y;
						nu[nb * 3 + 2] = 0;
						norm_v(&nu[nb * 3]);
						//infinite, define qtensor and don't evolve any more
						//homeotropic noninfinite, define qo
						if(infinite == 1){
							for(n = 0; n < 6; n ++){
								Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S);
							}		
						}
						else if(degenerate == 0 && infinite == 0){
							for(n = 0; n < 6; n ++){
								Qo[nb * 6 + n] = dir2ten(&nu[nb * 3], n, S);
					//			Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S);
							}		
						}
					//	else if(degenerate == 1){
					//		for(n = 0; n < 6; n ++){
					//			Qold[nd * 6 + n] = dir2ten(dir, n, S);
					//		}		
					//	}
					}	
					else if(nboundary[l]){
						for(m = 0; m < Np; m ++){
							x = (i- pos[m][0])*dx ;
							y = (j- pos[m][1])*dy ;
							z = (k- pos[m][2])*dz ;
							if(z > 0.5 * Nz){
								z -= Nz;
							}
							else if(z < -0.5 * Nz){
								z += Nz;
							}
							dis = x*x+y*y+z*z;
							if (dis <= (Rp + 0.5) * (Rp + 0.5)){
								if (dis == 0){
									printf("Error in neighbors on particle boundary.\n");
									return false;
								}
								else {
									nu[nb * 3 + 0] = x;
									nu[nb * 3 + 1] = y;
									nu[nb * 3 + 2] = z;
									norm_v(&nu[nb * 3]);
								}
								if(pos[m][3] == 0){
									//Infinite, do not evolve;
									share[nd] = 8;
									for(n = 0; n < 6; n ++){
										Qold[nd * 6 + n] = dir2ten(&nu[nb * 3], n, S);
									}		
								}
								else if(pos[m][3] == 1){
									//Non-Inf H
									share[nd] = 4;
									for(n = 0; n < 6; n ++){
										Qo[nb * 6 + n] = dir2ten(&nu[nb * 3], n, S);
									}		
								}
								else if(pos[m][3] == 2){
									//Degenerate
									share[nd] = 6;
								}
								break;
							}
						}
						nsurf --;
					}
					//define boundary
					if(nu[nb * 3 + 0] >= 0){
						neighbor[nd * 6 + 0] = indx[i + 1 + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = indx[i + 2 + j * Nx + k * Nx * Ny];
					}
					else if(nu[nb * 3 + 0] < 0){
						neighbor[nd * 6 + 0] = indx[i - 1 + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = indx[i - 2 + j * Nx + k * Nx * Ny];
					}
					if(nu[nb * 3 + 1] >= 0){
						neighbor[nd * 6 + 2] = indx[i + (j + 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = indx[i + (j + 2) * Nx + k * Nx * Ny];
					}
					else if(nu[nb * 3 + 1] < 0){
						neighbor[nd * 6 + 2] = indx[i + (j - 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = indx[i + (j - 2) * Nx + k * Nx * Ny];
					}
					if(nu[nb *3 + 2] >= 0){
						neighbor[nd * 6 + 4] = indx[i + j * Nx + peri(k + 1, 2) * Nx * Ny];
						neighbor[nd * 6 + 5] = indx[i + j * Nx + peri(k + 2, 2) * Nx * Ny];
					}
					else if(nu[nb * 3 + 2] < 0){
						neighbor[nd * 6 + 4] = indx[i + j * Nx + peri(k - 1, 2) * Nx * Ny];
						neighbor[nd * 6 + 5] = indx[i + j * Nx + peri(k - 2, 2) * Nx * Ny];
					}
					nb ++;
				}
				l ++;
			}
		}
	}

	if (nb != surf){
		printf("Problem in initialization of share. nb is %d not equal to surf %d.\n", nb, surf);
		return false;
	}

	for(nd = 0; nd < droplet; nd ++){
		//for all Bulk point, if one of the neighbor is surface point
		count1 = 0;
		if(share[nd] == 0){
			for(n = 0; n < 6; n ++){
				if(share[neighbor[nd * 6 + n]] >= 2){
					count1 ++;
				}
			}
			if(count1 > 1){
				share[nd] += 1;
			} 
		}
		//for all surface point, if one of the neighbor is not defined
		else if(share[nd] < 8 && share[nd] >= 2){	
			for(n = 0; n < 6; n++){
				if(neighbor[nd * 6 + n] == -1){
					count1 ++;	
				}
			}
			if(count1 > 0){
				share[nd] += 1;
			} 
		}
		//for all nodes with problem, share +1
	}			
	free(ndrop);
	free(indx);
	if(Np != 0){
		for(m = 0; m < Np; m ++)	free(pos[m]);
		free(pos);
	}
	printf("Initialization of droplet successfull.\n");
	return true;
}

