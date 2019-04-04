#ifndef PARALLELSETUP_H
#define PARALLELSETUP_H

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
	double * data1d;
	double ** data2d;

	double* halos1d;
	double**	halos2d;
	int id;
	int nbrup, nbrdown, nbrleft, nbrright;
	int size;

} Block;

void create_vec(Block *block, MPI_Datatype *vec);
void halo_swaps(Block* block, MPI_Comm comm, MPI_Datatype vec);
void find_my_neighbours(Block* my_block);
void free_block(Block* my_block);
void find_my_neighbours(Block* my_block);
void alloc_block(Block* my_block, int size);

#endif
