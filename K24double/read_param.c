#include "finite.h"

//read parameters from file param.in
bool read_param() {
        FILE *param;
	param = fopen("param.in","r");
	if(param == (FILE*)NULL){
		printf("File param.in not found.\n");
		return false;
	}
        fscanf(param,"Nx %d\n", &Nx);
        fscanf(param,"Ny %d\n", &Ny);
        fscanf(param,"Nz %d\n", &Nz);
        fscanf(param,"Lx %lf\n", &Lx);
        fscanf(param,"Ly %lf\n", &Ly);
        fscanf(param,"Lz %lf\n", &Lz);
        fscanf(param,"W %lf\n", &W);
        fscanf(param,"U %lf\n", &U);
        fscanf(param,"L1 %lf\n", &L1);
        fscanf(param,"L2 %lf\n", &L2);
        fscanf(param,"L3 %lf\n", &L3);
        fscanf(param,"L4 %lf\n", &L4);
        fscanf(param,"chiral %d\n", &chiral);
        fscanf(param,"N %lf\n", &N);
	//geo = 0 : bulk; geo = 1: channel; geo = 2: cylinder; geo = 3: drop; geo = 4: ellip
        fscanf(param,"geo %d\n", &geo);
        fscanf(param,"degenerate %d\n", &degenerate);
        fscanf(param,"infinite %d\n", &infinite);
        fscanf(param,"Np %d\n", &Np);
        fscanf(param,"Rp %lf\n", &Rp);
        fscanf(param,"Wp %lf\n", &Wp);
     //   fscanf(param,"pdegenerate %d\n", &pdegenerate);
       // fscanf(param,"pinfinite %d\n", &pinfinite);
        fscanf(param,"seed %d\n", &seed);
        fscanf(param,"rand_seed %d\n", &rand_seed);
        fscanf(param,"tmin tmax %lf %lf\n", &tmin, &tmax);
        fscanf(param,"increment %d\n", &increment);
        fscanf(param,"accuracy %lf\n", &accuracy);
        fscanf(param,"init_dir %lf %lf %lf\n", &init_dir[0], &init_dir[1], &init_dir[2]);
        fscanf(param,"dir1 %lf %lf %lf\n", &dir1[0], &dir1[1], &dir1[2]);
        fscanf(param,"dir2 %lf %lf %lf\n", &dir2[0], &dir2[1], &dir2[2]);
	if(myid == root){
		printf("Nx %d\n", Nx);
		printf("Ny %d\n", Ny);
		printf("Nz %d\n", Nz);
		printf("Lx %lf\n", Lx);
		printf("Ly %lf\n", Ly);
		printf("Lz %lf\n", Lz);
		printf("W %lf\n", W);
		printf("U %lf\n", U);
		printf("L1 %lf\n", L1);
		printf("L2 %lf\n", L2);
		printf("L3 %lf\n", L3);
		printf("L4 %lf\n", L4);
		printf("chiral %d\n", chiral);
		printf("N %lf\n", N);
		printf("geo %d\n", geo);
		printf("degenerate %d\n", degenerate);
		printf("infinite %d\n", infinite);
		printf("Np %d\n", Np);
		printf("Rp %lf\n", Rp);
		printf("Wp %lf\n", Wp);
	//	printf("pdegenerate %d\n", pdegenerate);
	//	printf("pinfinite %d\n", pinfinite);
		printf("seed %d\n", seed);
		printf("rand_seed %d\n", rand_seed);
		printf("tmin tmax %lf %lf\n", tmin, tmax);
		printf("increment %d\n", increment);
		printf("accuracy %lf\n", accuracy);
		printf("init_dir is %lf,%lf,%lf\n",init_dir[0],init_dir[1],init_dir[2]);
		printf("dir1 is %lf,%lf,%lf\n",dir1[0],dir1[1],dir1[2]);
		printf("dir2 is %lf,%lf,%lf\n",dir2[0],dir2[1],dir2[2]);
		if(infinite && degenerate){
			printf("Degenerate planar anchoring cannot be infinite.\n");    
			return false;
		}
		if(!norm_v(dir1) || !norm_v(dir2)){
			printf("Invalid dir1 or dir2 input.\n");
			return false;
		}	
		if(tmax < tmin){
			printf("Error for time input.\n");
			return false;
		}
	}
        fclose(param);
        return true;
}

bool read_nppos(double **pos){
	int Npnum, i, j;
	int coor;
	double x1, y1, z1, x2, y2, z2;
        FILE *nppos;
	nppos = fopen("nppos.in","r");
	if(nppos == (FILE*)NULL){
		printf("File nppos.in not found.\n");
		return false;
	}
	//read in NP number Np 
	//read in Coordinate: 0 for cartesian; 1 for spherical
        fscanf(nppos, "Np\t%d\nCoor\t%d\n", &Npnum, &coor);
	if(Npnum != Np){
		printf("Error in nanoparticle position file (NP number %d and %d conflict.)\n", Np, Npnum);
		return false;
	}
	if(coor != 0 && coor != 1){
		printf("Invalid input for nppos.in.\nSpherical coordinate: input 1; Cartesian coordinate: input 0.\n");
		return false;
	}
	fscanf(nppos, "posx\tposy\tposz\tanchoring\n");
	//read in particle position and anchoring(pos[i][3]), if spherical coordinate, change to cartesian
	//anchoring: 
	//0 for infinite anchoring
	//1 for noninfinite homeotropic
	//2 for planar degenerate
	for(i = 0; i < Np; i ++){
		fscanf(nppos, "%lf\t%lf\t%lf\t%lf\n", &pos[i][0], &pos[i][1], &pos[i][2], &pos[i][3]);
		if(pos[i][3] > 2 || pos[i][3] < 0){
			printf("Wrong input for anchoring condition.\n");
			return false;
		}
		if(coor == 1){
			x1 = lrint(Nx * 0.5) + pos[i][0] * sin(pos[i][1] / 180.0 * M_PI) * cos(pos[i][2] / 180.0 * M_PI);
			y1 = lrint(Ny * 0.5) + pos[i][0] * sin(pos[i][1] / 180.0 * M_PI) * sin(pos[i][2] / 180.0 * M_PI);
			z1 = lrint(Nz * 0.5) + pos[i][0] * cos(pos[i][1] / 180.0 * M_PI);
			pos[i][0] = x1;
			pos[i][1] = y1;
			pos[i][2] = z1;
		}
/*		if((pos[i][2] + Rp > Nz - 2) || (pos[i][2] - Rp < 2)){
			printf("Particle colloids with the wall.\n");
			fclose(nppos);
			return false;
		}
		if((pos[i][1] + Rp > Ny - 2) || (pos[i][1] - Rp < 2)){
			printf("Particle colloids with the wall.\n");
			fclose(nppos);
			return false;
		}
		if((pos[i][0] + Rp > Nx - 2) || (pos[i][0] - Rp < 2)){
			printf("Particle colloids with the wall.\n");
			fclose(nppos);
			return false;
		}
*/
		printf("%lf\t%lf\t%lf\t%lf\n", pos[i][0], pos[i][1], pos[i][2], pos[i][3]);
	}
	//examine if particles overlap
	if(Np > 1){
		for(i = 0; i < Np - 1; i ++){
			for(j = i + 1; j < Np; j ++){
				x1 = pos[i][0] - pos[j][0];
				y1 = pos[i][1] - pos[j][1];
				z1 = pos[i][2] - pos[j][2];
				if((x1 * x1) + (y1 * y1) + (z1 * z1) <= (4 * (Rp + 2) * (Rp + 2))){
					printf("Particle %d and %d overlap.\n", i, j);
					printf("Distance is %lf.\n\n", sqrt(x1 * x1 + y1 * y1 + z1 * z1));
					return false;
				}
			}
		}
	}	

	//print nanoparticle positions
	fclose(nppos);
	return true;
}
