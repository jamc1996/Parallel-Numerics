#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "hhorth.h"
#include "parallelsetup.h"
#include "dense.h"
#include "tsqr.h"
#include "caqr.h"

/*	main.c -- program to run/test the tsqr and caqr functions implemented else where
 *
 *	Author: John Cormican
 *
 *	Purpouse: testing and timing of functions
 *
 *	Usage: compile with other programs in repository and run to perform parallel qr factorization
 *
 * */


void run_test(int n, int m, int b, int nprocs, int myid, MPI_Comm comm, int nbrup, int nbrdown);


int main(int argc, char* argv[])
/* Main function for testing the efficiency of parallel qr decomposition. */
{
	int myid, nprocs;
	MPI_Comm cartcomm;												// Variables for the Cartesian topology
	int ndims, dims[1], periods[1], reorder;
	int nbrdown, nbrup;
	int b = 50;
	int n = 10000;
	int m = 5000;
	int c;
	while((c = getopt(argc,argv,"n:m:b:")) != -1)
	{
    switch(c)
		{
		case 'n':
			n = atoi(optarg);
			break;
		case 'm':
			m = atoi(optarg);
			break;
		case 'b':
			b = atoi(optarg);
			break;
		case '?':
			return 1;  	
		}
	}

	MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Values set for 2d Cartesian topology:
	ndims      = 1;
  dims[0]    = nprocs;
  periods[0] = 0;
  reorder    = 0;

	// cartesian communicator created:
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);
  MPI_Cart_shift(cartcomm, 0, 1, &nbrup, &nbrdown);

	run_test(n, m, b, nprocs, myid, cartcomm, nbrup, nbrdown);


	MPI_Finalize();

	return 0;
}


void run_test(int n, int m, int b, int nprocs, int myid, MPI_Comm comm, int nbrup, int nbrdown)
/*Function to run the communication avoiding tall skinny qr factorization*/
{
	clock_t start, end;
	double time;

	DenseMatrix A = CreateNullMatrix(n,m);
	DenseMatrix R = CreateNullMatrix(n,m);

	if (myid == 0){
		SetToRandom(&A);
	}

 	MPI_Bcast(A.data_, n*m, MPI_DOUBLE, 0, comm);

	memcpy(R.data_,A.data_,sizeof(double)*n*m);

	MPI_Barrier(comm);
	if (myid == 0)
	{
		start = clock();
	}

	caqr(&R, nprocs, myid, b, comm, nbrup, nbrdown);

	MPI_Barrier(comm);
	if (myid == 0)
	{
		end = clock();
		time = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\nOn %d processors, using block size %d, communication avoiding qr factorization of a %d x %d matrix, took %lf seconds.\n\n",nprocs, b, n, m, time);
	}


	FreeMatrix(&A);
	FreeMatrix(&R);

}
