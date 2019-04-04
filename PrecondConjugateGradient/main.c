#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "problemsetup.h"
#include "vectorops.h"
#include "conjugateg.h"
#include "parallelsetup.h"
#include "precondition.h"

int main(int argc, char *argv[])
{
	int block_size = 5;
	int myid, nprocs=1;
	//MPI_Init(&argc, &argv);
  //MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	//MPI_Comm_size(MPI_COMM_WORLD, &nprocs);		
	MPI_Datatype vec;

	if(nprocs == 1)
	{
		SparseMatrix L;
		Grid x, b;
		allocate_L(&L, &x);
		alloc_serial_grid(&x, block_size);
		alloc_serial_grid(&b, block_size);
		fill_grid(&x);

		define_b(&x, &b);

		//rb_gauss_seidel(&x, &b);

		conjugate_gradient(&x, &b);

		//allocate_L(&L, &x);

		print_full_grid(&x);
	
		free_grid(&x);
	}
	if(nprocs == 8)
	{
		Block my_block;
		my_block.id = myid;
		find_my_neighbours(&my_block);
		alloc_block(&my_block,block_size);
		create_vec(&my_block, &vec);
		//halo_swaps(&my_block, MPI_COMM_WORLD, vec);
		printf("myid is %d and my neighbours are %d %d %d %d\n",myid,my_block.nbrup,my_block.nbrdown,my_block.nbrleft,my_block.nbrright);

		free_block(&my_block);

	}
	//MPI_Finalize();	
	return 0;
}














