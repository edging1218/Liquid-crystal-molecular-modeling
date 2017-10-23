#include "finite.h"

int peri(int node, int dir){
	/*
	Periodic boundary condition
	*/
	if(dir == 0){
		if(node >= 0 && node < Nx)	
			return node;
		else if(node >= Nx)
			return node - Nx;
		else
			return node + Nx;
	}
	else if(dir == 1){
		if(node >= 0 && node < Ny)	
			return node;
		else if(node >= Ny)
			return node - Ny;
		else 
			return node + Ny;
	}
	else if(dir == 2){
		if(node >= 0 && node < Nz)	
			return node;
		else if(node >= Nz)
			return node - Nz;
		else 
			return node + Nz;
	}
	else{
		printf("Wrong input for function periodic.\n");
	}
	return 0;
}

bool initial(){
	/*
	Initialize the system with given geometry. 
	Write the header of energy output file.
	Write the grid info to grid.bin file.
	*/
	int signal, l;
	//geo = 0: bulk
	//geo = 1: channel
	//geo = 2: drop
	//geo = 3: ellip
	if(geo == 0){
		if(!initial_bulk()){
			return false;
		}
	}
	else if(geo == 1){
		if(!initial_channel()){
			return false;
		}
	}
	else if(geo == 2){
		if(!initial_cylinder()){
			return false;
		}
	}
	else if(geo == 3){
		if(!initial_drop()){
			return false;
		}
	}
	else if(geo == 4){
		if(!initial_ellip()){
			return false;
		}
	}

	printf("\nS is %lf.\n\n", S);

	FILE* energy;
	FILE* grid;
	energy = fopen("energy.out", "w");
	fprintf(energy,"cycle\tEnergy_diff\tEnergy_ldg\tEnergy_l1\tEnergy_l2\tEnergy_l3\tEnergy_l4\tEnergy_chiral\tEnergy_surf\tEnergy_surf\tEnergy_tot\n");
	fclose(energy);

	grid = fopen("grid.bin", "wb");
	l = 0;
	for(l = 0; l < tot; l++){
		if(boundary[l] || nboundary[l])	signal = 1;
		else if(drop[l]) 	signal = 0;			
		else signal = -1;
		fwrite(&signal, sizeof(int), 1, grid);
	}
	fclose(grid);	
	return true;
}
