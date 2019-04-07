#ifndef PRECONDITION_H
#define PRECONDITION_H

#include "structs.h"

#include "problemsetup.h"
#include "parallelsetup.h"
#include "haloswaps.h"


void rb_gauss_seidel(Grid* p, Grid* r);
void par_rb_gauss_seidel(Block* p, Block* r, MPI_Comm comm, MPI_Datatype vec);

#endif
