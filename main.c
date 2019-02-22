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
{
	int myid, nprocs;

	int b = 4;
	int n = 200;
	int i,j;
	int m = 20;

	DenseMatrix A = CreateNullMatrix(n,b);
	DenseMatrix Rf = CreateNullMatrix(n,b);
	DenseMatrix Qf = CreateNullMatrix(n,b); 	

	DenseMatrix A2 = CreateNullMatrix(n,m);
	DenseMatrix Rf2 = CreateNullMatrix(n,m);
	DenseMatrix Qf2 = CreateNullMatrix(n,m); 	

	MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (myid == 0){
		SetToRandom(&A);
		SetToRandom(&A2);
	}
	
 	MPI_Bcast(A.data_, n*b, MPI_DOUBLE, 0,MPI_COMM_WORLD);
 	MPI_Bcast(A2.data_, n*m, MPI_DOUBLE, 0,MPI_COMM_WORLD);
	//tsqr(&A, &Rf, &Qf, nprocs, myid, b);
	memcpy(Rf2.data_,A2.data_,sizeof(double)*n*m);
		

/*	if(myid == 0){printf("Rr[0][%d] = %lf\n",0,Rf2.entry[0][0]);
		for (i=0;i< 15;i++){
			for (j=0; j< 15; j++){
				printf("  %lf  ",A2.entry[j][i]);
			} printf("\n");}
	}*/

	caqr(&Rf2, &Qf2, nprocs, myid, b);
	
	//Find Qt A
	

	if (myid == 2){
		DenseMatrix Rr = CreateNullMatrix(n,m);
		DenseMatrix V = CreateNullMatrix(n,m);
		hhorth(&A2,&V,&Rr);
		sleep(1);

		printf("Rr[0][%d] = %lf\n",0,Rr.entry[0][0]);
		for (i=0;i< 15;i++){
			for (j=0; j< 15; j++){
				printf("  %lf  ",Rr.entry[j][i]);
			}printf("\n");}
	}
	FreeMatrix(&A2);

	MPI_Finalize();

	return 0;
}
