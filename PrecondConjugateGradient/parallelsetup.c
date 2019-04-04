#include "parallelsetup.h"

void alloc_block(Block* my_block, int size)
{
	my_block->data1d = malloc(sizeof(double)*size*size);
	my_block->data2d = malloc(sizeof(double*)*size);

	my_block->halos1d = malloc(sizeof(double)*size*4);
	my_block->halos2d = malloc(sizeof(double*)*4);

	int i;
	for ( i=0; i<size; i++)
	{
		my_block->data2d[i] = &my_block->data1d[i*size];
	}

	for ( i=0; i<4; i++)
	{
		my_block->halos2d[i] = &my_block->halos1d[i*size];
	}
}

void create_vec(Block *block, MPI_Datatype *vec)
{
	MPI_Type_vector(block->size, 1, block->size, MPI_DOUBLE, vec);
	MPI_Type_commit(vec);
}


void halo_swaps(Block* block, MPI_Comm comm, MPI_Datatype vec)
{
	if(block->nbrup > 0 && block->nbrdown > 0)
	{
		// Send up, receive from below
		MPI_Sendrecv(block->data2d[0], block->size, MPI_DOUBLE, block->nbrup, 0, block->halos2d[2], block->size, MPI_DOUBLE, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);

		// Send down, receive from above
		MPI_Sendrecv(block->data2d[block->size-1], block->size, MPI_DOUBLE, block->nbrdown, 0, block->halos2d[0], block->size, MPI_DOUBLE, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
	}
	else if(block->nbrup > 0)
	{
		// Send UP
		MPI_Send(block->halos2d[0], block->size, MPI_DOUBLE, block->nbrup, 0, comm);

		// Receive from above
		MPI_Recv(block->data2d[0], block->size, MPI_DOUBLE, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
	}
	else if(block->nbrdown > 0)
	{
		// Receive from below
		MPI_Recv(block->data2d[block->size-1], block->size, MPI_DOUBLE, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);
		// Send down
		MPI_Send(block->halos2d[2], block->size, MPI_DOUBLE, block->nbrdown, 0, comm);
	}



}

void free_block(Block* my_block)
{
	free(my_block->data1d);
	free(my_block->data2d);

	free(my_block->halos1d);
	free(my_block->halos2d); 
}

void find_my_neighbours(Block* my_block)
{
	if(my_block->id == 0)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = 3;	
		my_block->nbrright = 1;
		my_block->nbrleft = -1;
	}
	if(my_block->id == 1)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = -1;	
		my_block->nbrright = 2;
		my_block->nbrleft = 1;
	}	
	if(my_block->id == 2)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = -1;	
		my_block->nbrright = -2;
		my_block->nbrleft = 1;
	}	
	if(my_block->id == 3)
	{
		my_block->nbrup = 0;
		my_block->nbrdown = 4;	
		my_block->nbrright = -1;
		my_block->nbrleft = -1;
	}	
	if(my_block->id == 4)
	{
		my_block->nbrup = 3;
		my_block->nbrdown = 6;	
		my_block->nbrright = 5;
		my_block->nbrleft = -1;
	}	
	if(my_block->id == 5)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = 7;	
		my_block->nbrright = -1;
		my_block->nbrleft = 4;
	}	
	if(my_block->id == 6)
	{
		my_block->nbrup = 4;
		my_block->nbrdown = -3;	
		my_block->nbrright = 7;
		my_block->nbrleft = -1;
	}	
	if(my_block->id == 7)
	{
		my_block->nbrup = 5;
		my_block->nbrdown = -3;	
		my_block->nbrright = -1;
		my_block->nbrleft = 6;
	}	
}
