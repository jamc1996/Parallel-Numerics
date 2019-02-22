#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#include "hhorth.h"
#include "parallelsetup.h"
#include "dense.h"
#include "tsqr.h"
#include "caqr.h"

int main(int argc, char* argv[])
/* Main function for testing the efficiency of parallel qr decomposition. */
{
	int myid, nprocs;

	int b = 4;
	int n = 200;
	int m = 20;

	int i,j;

	DenseMatrix Ab = CreateNullMatrix(n,b);
	DenseMatrix Rb = CreateNullMatrix(n,b);
	DenseMatrix Qb = CreateNullMatrix(n,b); 	

	DenseMatrix Am = CreateNullMatrix(n,m);
	DenseMatrix Rm = CreateNullMatrix(n,m);
	DenseMatrix Qm = CreateNullMatrix(n,m); 	

	MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (myid == 0){
		SetToRandom(&Ab);
		SetToRandom(&Am);
	}
	
 	MPI_Bcast(Ab.data_, n*b, MPI_DOUBLE, 0,MPI_COMM_WORLD);
 	MPI_Bcast(Am.data_, n*m, MPI_DOUBLE, 0,MPI_COMM_WORLD);
	
	memcpy(Rb.data_,Ab.data_,sizeof(double)*n*b);
	tsqr(&Rb, &Qb, nprocs, myid, b);

	//memcpy(Rm.data_,Am.data_,sizeof(double)*n*m);
	//caqr(&Rm, &Qm, nprocs, myid, b);
	
	//Find Qt A
	

/*	if (myid == 2){
		hhorth(&A2,&V,&Rr);
		sleep(1);

		printf("Rr[0][%d] = %lf\n",0,Rr.entry[0][0]);
		for (i=0;i< 15;i++){
			for (j=0; j< 15; j++){
				printf("  %lf  ",Rr.entry[j][i]);
			}printf("\n");}
	}*/

	FreeMatrix(&Ab);
	FreeMatrix(&Rb);
	FreeMatrix(&Qb);
	FreeMatrix(&Am);
	FreeMatrix(&Rm);
	FreeMatrix(&Qm);

	MPI_Finalize();

	return 0;
}
