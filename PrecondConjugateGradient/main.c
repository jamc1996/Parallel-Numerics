#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

#include "problemsetup.h"
#include "vectorops.h"
#include "conjugateg.h"
#include "parallelsetup.h"
#include "precondition.h"
#include "haloswaps.h"

/*	main.c -- program for running and testing code for 
 *	the conjugate gradient method in serial and parallel.
 *
 *	Author: John Cormican
 *
 *	Purpouse: Testing of implementation of preconditioned conjugate gradient 
 *
 *	Usage: To be run either basically or with mpi exec -n 8
 */


int main(int argc, char *argv[])
/* Function to test code for performance of parallel methods. */
{
	int c;	
	int print_flag = 0;
	int block_size = 4;
	while ((c = getopt(argc, argv, "pb:")) != -1)
	{
		switch (c)
		{
		case 'p':
			print_flag = 1;
			break;
		case 'b':
			block_size = atoi(optarg);
			break;
		case '?':
			return 1;
		}
	}

	int myid, nprocs=1;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);		
	MPI_Datatype vec;
	
	// Serial implementation run:
	if(nprocs == 1)
	{
		// x and b grids allocated and filled for mat eqn Ax=b
		Grid x, b;
	
		alloc_serial_grid(&x, block_size);
		alloc_serial_grid(&b, block_size);

		fill_grid(&x);
		define_b(&x, &b);

		// Preconditioned conjugate gradient used to solve:
		precond_conjugate_gradient(&x, &b);

		if(print_flag == 1)
		{
			print_full_grid(&x);
		}
		free_grid(&x);
		free_grid(&b);
	}
	// Parallel implementation run:
	else if(nprocs == 8)
	{
		// x and b grids allocated and filled for mat eqn Ax=b
		Block x, b;

		alloc_block(&x,block_size, myid);
		alloc_block(&b,block_size, myid);

		initiate_block(&x,MPI_COMM_WORLD);
		fill_blockB(&b);

		// MPI datatype for passing non-contiguous data
		create_vec(&x, &vec);	

		parallel_precond_cg(&x, &b, MPI_COMM_WORLD, vec);
				
		MPI_Barrier(MPI_COMM_WORLD);
		if(myid == 0)
		{
			printf("x is now\n");
		}
		print_block(&x,MPI_COMM_WORLD);

		free_block(&x);
		free_block(&b);
	}
	else
	{
		printf("main.c: please use either 1 or 8 processors.\n");
		MPI_Finalize();	
		return 1;
	}

	MPI_Finalize();	
	return 0;
}














