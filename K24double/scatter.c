#include "finite.h"

bool scatter(){
	int i, j;
	int count;
	int count_tot;
	int* count_root;
	int* displ;

	MPI_Bcast(&idx, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&idy, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&idz, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&iddx, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&iddy, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&iddz, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&qch, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&dV, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&dAdrop, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&dApart, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(&droplet, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(&length, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(&AnchNInf, 1, MPI_BYTE, root, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	//define shared window and store Qold of root processor to q for all processors to access
        MPI_Win_allocate_shared(6 * length * sizeof(double), 1, MPI_INFO_NULL, shmcomm, &q, &win);
        MPI_Scatter(Qold, 6 * length, MPI_DOUBLE, q, 6 * length, MPI_DOUBLE, root, MPI_COMM_WORLD);

	//define shared window and store neighbor of root processor to neigb for all processors to access
        MPI_Win_allocate_shared(6 * length * sizeof(int), 1, MPI_INFO_NULL, shmcomm, &neigb, &win2);
        MPI_Scatter(neighbor, 6 * length, MPI_INT, neigb, 6 * length, MPI_INT, root, MPI_COMM_WORLD);

	//populate share at root processor to sign at all processors
	sign = (int*)malloc(length * sizeof(int));
	for(i = 0; i < length; i ++)	sign[i] = -1;
        MPI_Scatter(share, length, MPI_INT, sign, length, MPI_INT, root, MPI_COMM_WORLD);
	
	//Allocate Qnew(qn)
	qn = (double*)malloc(6 * length * sizeof(double));
	for(i = 0; i < 6 * length; i ++)	qn[i] = q[i];	

	//Verify the number of droplet and boundary. If not consistent, report error.
	count = 0;
	for(i = 0; i < length; i ++){
		if(sign[i] >= 0 && sign[i] < 10)	count ++;
	}

 	MPI_Reduce(&count, &count_tot, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
	if(myid == root && count_tot != droplet){
		printf("Error in scatter. Counted number %d is not equal to droplet %d.\n", count_tot, droplet);	
		return false;
	}
	count = 0;
	for(i = 0; i < length; i ++){
		if(sign[i] >= 2 && sign[i] < 10)	count ++;
	}
 	MPI_Reduce(&count, &count_tot, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
	if(myid == root && count_tot != surf){
		printf("Error in scatter(boundary). Counted number %d is not equal to surface %d.\n", count_tot, surf);	
		return false;
	}
	count *= 3;

	nu_p = (double*)malloc(count * sizeof(double));
	count_root = (int*)malloc(numprocs * sizeof(int));
	displ = (int*)malloc(numprocs * sizeof(int));

//	if(myid == root) printf("Check3.\n");
//	scatter nu and qo defined at boundary nodes to different processors.
	MPI_Gather(&count, 1, MPI_INT, count_root, 1, MPI_INT, root, MPI_COMM_WORLD);		
	if(myid == root){
		for(i = 0; i < numprocs; i ++){
			displ[i] = 0;
			for(j = 0; j < i; j++){
				displ[i] += count_root[j];
			}
		}
	}
        MPI_Scatterv(nu, count_root, displ, MPI_DOUBLE, nu_p, count, MPI_DOUBLE, root, MPI_COMM_WORLD);

	
	if((degenerate == 0 && infinite == 0) || AnchNInf){
		count *= 2;
		if(myid == root){	
			for(i = 0; i < numprocs; i ++){
				count_root[i] *= 2;
				displ[i] *= 2;
			}
		}
		qo_p = (double*)malloc(count * sizeof(double));
		MPI_Scatterv(Qo, count_root, displ, MPI_DOUBLE, qo_p, count, MPI_DOUBLE, root, MPI_COMM_WORLD);
	}

//	printf("check4.\n");
	if(myid == root){
		free(neighbor);
		free(Qold);
		free(share);
		free(nu);
		if((degenerate == 0 && infinite == 0) || AnchNInf)	free(Qo);
		printf("Scattering successful.\n");
	}
	free(count_root);
	free(displ);
	return true;

}





