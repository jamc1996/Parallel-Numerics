#ifndef HALOSWAPS_H
#define HALOSWAPS_H

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

#include "structs.h"


void create_vecs(Block *block, MPI_Datatype *vecA, MPI_Datatype *vecB, MPI_Datatype *vecC, MPI_Datatype *vecD);
void create_vec(Block *block, MPI_Datatype *vec);
void create_rows(Block *block, MPI_Datatype *vecA, MPI_Datatype *vecB);
void halo_swaps(Block* block, MPI_Comm comm, MPI_Datatype vecA, MPI_Datatype vecB, MPI_Datatype vecC, MPI_Datatype vecD, MPI_Datatype rowA, MPI_Datatype rowB, int oddr, int oddv);
void halo_swapping(Block* block, MPI_Comm comm, MPI_Datatype vec);
#endif
