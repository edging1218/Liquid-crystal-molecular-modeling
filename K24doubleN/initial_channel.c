#include "finite.h"

bool initial_channel(){
	int l;
	int nb, nd, np;
	int i, j, k, n, m;
	double dis, x, y, z;
	int xm, xp, ym, yp, zm, zp;
	int nsurf;
	bool *ndrop;
	int *indx;
	double **pos;
	int count1;
	double dir[3] = {1, 0, 0};

	double dx = Lx/(Nx-1);
	double dy = Ly/(Ny-1);
	double dz = Lz/(Nz-1);

	idx = 1 / dx;
	idy = 1 / dy;
	idz = 1 / dz;
	iddx = idx * idx;
	iddy = idy * idy;
	iddz = idz * idz;

        qch = 2 * M_PI * N / Nz;

	surf = 2 * Nx * Ny;
	nsurf = 0;
	tot = Nx * Ny * Nz;
	bulk = Nx * Ny * (Nz - 2);

	//allocate drop and boundary
	ndrop = (bool*)malloc(tot * sizeof(bool));
	nboundary = (bool*)malloc(tot * sizeof(bool));
	drop = (bool*)malloc(tot * sizeof(bool));
	boundary = (bool*)malloc(tot * sizeof(bool));
	indx = (int*)malloc(tot * sizeof(int));
	for(l = 0; l < tot; l ++){
		drop[l] = true;
		boundary[l] = false;
		ndrop[l] = false;
		nboundary[l] = false;
		indx[l] = -1;
	}

	//define the channel surface 
	for (j = 0; j < Ny; j++){
		for (i = 0; i < Nx; i++){
			for(k = 0; k <= Nz - 1; k += Nz - 1){
				l = i + j * Nx + k * Nx * Ny;
				drop[l] = false;
				boundary[l] = true;
			}
		}
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
					for(m = 0; m < Np; m ++){
						x = i - pos[m][0];
						y = j - pos[m][1];
						z = k - pos[m][2];
						if(x > 0.5 * Nx){
							x -= Nx;
						}
						else if(x < -0.5 * Nx){
							x += Nx;
						}
						if(y > 0.5 * Ny){
							y -= Ny;
						}
						else if(y < -0.5 * Ny){
							y += Ny;
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
							break;
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
					if(drop[l]){
						xm = peri(i - 1, 0) + j * Nx + k * Nx * Ny;
						xp = peri(i + 1, 0) + j * Nx + k * Nx * Ny;
						ym = i + peri(j - 1, 1) * Nx + k * Nx * Ny;
						yp = i + peri(j + 1, 1) * Nx + k * Nx * Ny;
						zm = i + j * Nx + (k - 1) * Nx * Ny;
						zp = i + j * Nx + (k + 1) * Nx * Ny;
						if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zm] || ndrop[zp]){
							nboundary[l] = true;
							drop[l] = false;
							bulk --;
							nsurf ++;	
							surf ++;
						}
					}
					else if(boundary[l] && Wp > W){
						xm = peri(i - 1, 0) + j * Nx + k * Nx * Ny;
						xp = peri(i + 1, 0) + j * Nx + k * Nx * Ny;
						ym = i + peri(j - 1, 1) * Nx + k * Nx * Ny;
						yp = i + peri(j + 1, 1) * Nx + k * Nx * Ny;
						if(k == 0)
							zp = i + j * Nx + (k + 1) * Nx * Ny;
						else
							zp = i + j * Nx + (k - 1) * Nx * Ny;
						if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zp]){
							boundary[l] = false;
							nboundary[l] = true;
							nsurf ++;
						}
					}
					l ++;
				}
			}
		}
	}
	//For particle on the surface
	//dV = (Nx * Ny * Nz - 4.0 / 3 * M_PI * Rp * Rp * Rp * Np * 0.5) / bulk; 
	//dAdrop = (2 * Nx * Ny - M_PI * Rp * Rp * Np) / (surf -  nsurf);
	//dApart = (2 * M_PI * Rp * Rp * Np) / nsurf ;
	
	//For particle in the bulk
	dV = (Lx * Ly * Lz - 4.0 / 3 * M_PI * Rp * Rp * Rp * Np) / bulk; 
	dAdrop = (2 * Lx * Ly) / (surf -  nsurf);
	if(Np != 0){
		dApart = (4 * M_PI * Rp * Rp * Np) / nsurf ;
	}
	else{
		dApart = 0;
	}
	printf("\ndV is %lf\ndA of droplet is %lf\ndA of nanoparticle is %lf\n", dV, dAdrop, dApart); 

	droplet = bulk + surf;
	printf("\nDroplet nodes number is %d.\nBulk nodes number is %d.\nDroplet surface nodes number is %d. \nParticle surface nodes number is %d.\n", droplet, bulk, surf, nsurf); 

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
		if(!ndrop[l]){
			indx[l] = nd;
			//bulk :  share/sign = 0
			//channel boundary :  share/sign = 2
			//NP boundary :  share/sign = 4
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
					neighbor[nd * 6 + 0] = indx[peri(i - 1, 0) + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 1] = indx[peri(i + 1, 0) + j * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 2] = indx[i + peri(j - 1, 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 3] = indx[i + peri(j + 1, 1) * Nx + k * Nx * Ny];
					neighbor[nd * 6 + 4] = indx[i + j * Nx + (k - 1) * Nx * Ny];
					neighbor[nd * 6 + 5] = indx[i + j * Nx + (k + 1) * Nx * Ny];
				}
				else if(nboundary[l] || boundary[l]){
					if(boundary[l]){
						if(k == 0){
							nu[nb * 3 + 0] = 0;
							nu[nb * 3 + 1] = 0;
							nu[nb * 3 + 2] = 1;
						}
						else if(k == Nz - 1){
							nu[nb * 3 + 0] = 0;
							nu[nb * 3 + 1] = 0;
							nu[nb * 3 + 2] = -1;
						}
						else{
							printf("Error in initializing neighbors for channel surface.\n");
							return false;
						}
						//infinite, define qtensor and don't evolve any more
						//homeotropic noninfinite, define qo
						if(infinite == 1){
							if(k == 0){
								for(n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(dir1, n, S);
								}	
							}	
							else{
								for(n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(dir2, n, S);
								}	
							}
						}
						else if(degenerate == 0 && infinite == 0){
							if(k == 0){
								for(n = 0; n < 6; n ++){
									Qo[nb * 6 + n] = dir2ten(dir1, n, S);
							//		Qold[nd * 6 + n] = dir2ten(dir1, n, S);
								}		
							}
							else{
								for(n = 0; n < 6; n ++){
									Qo[nb * 6 + n] = dir2ten(dir2, n, S);
							//		Qold[nd * 6 + n] = dir2ten(dir2, n, S);
								}		
							}
						}
/*
						else if(degenerate = 1){
							if(k == 0){
								for(n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}	
							}	
							else{
								for(n = 0; n < 6; n ++){
									Qold[nd * 6 + n] = dir2ten(dir, n, S);
								}	
							}
						}
*/
					}
					else{
						for(m = 0; m < Np; m ++){
							x = i - pos[m][0];
							y = j - pos[m][1];
							z = k - pos[m][2];
							if(x > 0.5 * Nx)
								x -= Nx;
							else if(x < -0.5 * Nx)
								x += Nx;
							if(y > 0.5 * Ny)
								y -= Ny;
							else if(y < -0.5 * Ny)
								y += Ny;
							dis = x*x+y*y+z*z;
							if (dis <= (Rp + 1) * (Rp + 1)){
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
						neighbor[nd * 6 + 0] = indx[peri(i + 1, 0) + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = indx[peri(i + 2, 0) + j * Nx + k * Nx * Ny];
					}
					else if(nu[nb * 3 + 0] < 0){
						neighbor[nd * 6 + 0] = indx[peri(i - 1, 0) + j * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 1] = indx[peri(i - 2, 0) + j * Nx + k * Nx * Ny];
					}
					if(nu[nb * 3 + 1] >= 0){
						neighbor[nd * 6 + 2] = indx[i + peri(j + 1, 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = indx[i + peri(j + 2, 1) * Nx + k * Nx * Ny];
					}
					else if(nu[nb * 3 + 1] < 0){
						neighbor[nd * 6 + 2] = indx[i + peri(j - 1, 1) * Nx + k * Nx * Ny];
						neighbor[nd * 6 + 3] = indx[i + peri(j - 2, 1) * Nx + k * Nx * Ny];
					}
					if(nu[nb *3 + 2] >= 0){
						neighbor[nd * 6 + 4] = indx[i + j * Nx + (k + 1) * Nx * Ny];
						neighbor[nd * 6 + 5] = indx[i + j * Nx + (k + 2) * Nx * Ny];
					}
					else if(nu[nb * 3 + 2] < 0){
						neighbor[nd * 6 + 4] = indx[i + j * Nx + (k - 1) * Nx * Ny];
						neighbor[nd * 6 + 5] = indx[i + j * Nx + (k - 2) * Nx * Ny];
					}
					nb ++;
				}
				l ++;
			}
		}
	}

	if(nsurf != 0){
		printf("Problems in initialization: nsurf is %d.\n", nsurf);
		return false;
	}

	for(nd = 0; nd < droplet; nd ++){
		//for andnd Bundk point, if one of the neighbor is surface point
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
	printf("Initialization of channel successfull.\n");
	return true;
}
