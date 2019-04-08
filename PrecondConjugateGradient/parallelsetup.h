#ifndef PARALLELSETUP_H
#define PARALLELSETUP_H

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

#include "problemsetup.h"
#include "structs.h"

void find_my_neighbours(Block* my_block);
void free_block(Block* my_block);
void find_my_neighbours(Block* my_block);
void alloc_block(Block* my_block, int size, int myid);
void initiate_block(Block* block, MPI_Comm comm);
void print_block(Block *block, MPI_Comm comm);
void fill_blockB(Block* blocked_b);
double pa_nbr_up(Block *b, int x, int y);
double pa_nbr_down(Block *b, int x, int y);
double pa_nbr_left(Block *b, int x, int y);
double pa_nbr_right(Block *b, int x, int y);

#endif
