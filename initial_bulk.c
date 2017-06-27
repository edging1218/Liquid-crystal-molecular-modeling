#include "finite.h"

bool initial_bulk(){
	int l;
	int nb, nd, np;
	int i, j, k, n, m;
	double dis, x, y, z;
	int xm, xp, ym, yp, zm, zp;
	bool *ndrop;
	int *indx;
	double **pos;
	int count1;

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

	tot = Nx * Ny * Nz;
	bulk = tot;

	//allocate drop and boundary
	ndrop = (bool*)malloc(tot * sizeof(bool));
	boundary = (bool*)malloc(tot * sizeof(bool));
	nboundary = (bool*)malloc(tot * sizeof(bool));
	drop = (bool*)malloc(tot * sizeof(bool));
	indx = (int*)malloc(tot * sizeof(int));
	for(l = 0; l < tot; l ++){
		drop[l] = true;
		boundary[l] = false;
		ndrop[l] = false;
		nboundary[l] = false;
		indx[l] = -1;
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
						if(z > 0.5 * Nz){
							z -= Nz;
						}
						else if(z < -0.5 * Nz){
							z += Nz;
						}
						dis = x*x+y*y+z*z;
						if (dis <= (Rp - 0.5) * (Rp - 0.5)){
							ndrop[l] = true;
							drop[l] = false;
							bulk --;
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
						zm = i + j * Nx + peri(k - 1, 2) * Nx * Ny;
						zp = i + j * Nx + peri(k + 1, 2) * Nx * Ny;
						if(ndrop[xm] || ndrop[xp] || ndrop[ym] || ndrop[yp] || ndrop[zm] || ndrop[zp]){
							nboundary[l] = true;
							drop[l] = false;
							bulk --;
							surf ++;
						}
					}
					l ++;
				}
			}
		}
	}
	//For particle in the bulk
	dV = (Lx * Ly * Lz - 4.0 / 3 * M_PI * Rp * Rp * Rp * Np) / bulk; 
	dAdrop = 0;
	if(Np > 0){
		dApart = (4 * M_PI * Rp * Rp * Np) / surf ;
	}
	else{
		dApart = 0;
	}
	printf("\ndV is %lf\ndA of nanoparticle is %lf\n", dV, dApart); 

	droplet = bulk + surf;
	printf("\nDroplet nodes number is %d.\nBulk nodes number is %d.\nSurface nodes number is %d.\n", droplet, bulk, surf); 

	//allocate nu 
	//allocate Qo only for finite homeotropic
	if(Np > 0){
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
		if(AnchNInf){
			Qo = (double*)malloc(6 * surf * sizeof(double));
			for(i = 0; i < surf * 6; i ++){
				Qo[i] = 0;
			}
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
			if(drop[l]) share[nd] = 0;
			else if(nboundary[l]){
				share[nd] = 4;
				nb ++;
			}
			else{
				printf("Error in initalization of grid.\n");
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
					neighbor[nd * 6 + 4] = indx[i + j * Nx + peri(k - 1, 2) * Nx * Ny];
					neighbor[nd * 6 + 5] = indx[i + j * Nx + peri(k + 1, 2) * Nx * Ny];
				}
				else if(nboundary[l]){
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
						if(z > 0.5 * Nz){
							z -= Nz;
						}
						else if(z < -0.5 * Nz){
							z += Nz;
						}
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

	for(l = 0; l < droplet; l++){
		count1 = 0;
		if(share[l] == 0){
			for(n = 0; n < 6; n ++){
				if(share[neighbor[l * 6 + n]] >= 2){
					count1 ++;
				}
			}
		}
		else if(share[l] >= 4 && share[l] < 8){	
			for(n = 0; n < 6; n++){
				if(neighbor[l * 6 + n] == -1){
					count1 ++;	
				}
			}
		}
		if(count1 != 0){
			share[l] += 1;
		} 
	}			

	free(ndrop);
	free(indx);
	if(Np != 0){
		for(m = 0; m < Np; m ++)	free(pos[m]);
		free(pos);
	}
	printf("Initalization of bulk successful.\n");
	return true;
}
