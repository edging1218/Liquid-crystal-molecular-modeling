#include <time.h>
#include "finite.h"

int main(int argc, char *argv[]){
        MPI_Init(&argc, &argv);
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	int i;
	bool flag = true;
	double deltat;

	double time_spend;
	double begin, end;

	FILE* energy;

	begin = MPI_Wtime();
	//read in parameters
	if(!read_param()){
		MPI_Comm_free(&shmcomm);
		MPI_Finalize();
		return 1;
	}
	
	dt = tmin;
	dE = 1;
	el_old = 1;
	cycle = 0;
	deltat = (tmax - tmin) / increment;

	S = 0.25 * (1 + 3 * sqrt(1 - 8 / (3 * U)));
	//	printf("Theoretical value is %lf.\n", third * (1 - third * U) * S * S - 2 * third * third * third * S * S * S * U + U / 9 * S * S * S * S );

	if(myid == root){
		//define droplet and boundary, introduce particles, initialize qtensor;
		if(!initial()){
			flag = false;
		}		
	}
	MPI_Bcast(&flag, 1, MPI_BYTE, root, MPI_COMM_WORLD);	
	
	//For all infomation initialized in root processor, scatter them to all other processors.
	if(!scatter()){
		flag = false;
	}		

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Win_fence(0, win);
	//Evolution
	while(flag){
		//Every 1000 steps calculate energies.
		if(cycle % 10000 == 0){
			free_energy();
			if(fabs(dE) < accuracy){
				flag = false;
			}
		}

		if(cycle % 10000 == 0){ 
			//Every 10000 steps check the trace of Qtensor
			if(myid == root){	
				for(i = 0; i < droplet; i++){
					//				checktr(&q[i * 6]);
					if(!checktr(&q[i * 6])){
						flag = false;
						printf("%d\n", i);
					}
				}
				if(!flag){
					printf("Error in the trace of q; cycle : %d.\n", cycle);
				}
				output();
			}	
			MPI_Bcast(&flag, 1, MPI_BYTE, root, MPI_COMM_WORLD);	
		}

		//Wait until all the processors are ready and relax the system, first bulk and then boundary.
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Win_fence(0, win);
		if(flag) relax_bulk();

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Win_fence(0, win);
		if(flag)	relax_surf();

		// Update learning rate
		if(dt < tmax){
			dt += deltat;
			if(dt >= tmax)	dt = tmax;
		}

		cycle ++;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Win_fence(0, win);
		
	}
	//Calculate free energy
	free_energy();

	end = MPI_Wtime();

	if(myid == root){
		//Write output file
		output();
		//Calculate time used.
		time_spend = (double)(end - begin) / 60.0;
		energy = fopen("energy.out", "a");
		if(time_spend < 60){
			fprintf(energy, "\nTime used:	%lf min.\n", time_spend);
			printf("\nTime used:	%lf min.\n", time_spend);
		}
		else{
			fprintf(energy, "\nTime used:	%lf h.\n", time_spend / 60.0);
			printf("\nTime used:	%lf h.\n", time_spend / 60.0);
		}
		fclose(energy);	
	}

	//deallocate dynamic arrays
	free_q();
        MPI_Win_free(&win);
        MPI_Win_free(&win2);
        MPI_Comm_free(&shmcomm);
        MPI_Finalize();

	return 0;
}
