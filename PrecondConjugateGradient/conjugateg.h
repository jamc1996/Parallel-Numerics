#ifndef CONJUGATEG_H
#define CONJUGATEG_H

#include "problemsetup.h"
#include "vectorops.h"
#include "precondition.h"
#include "parallelsetup.h"
#include "structs.h"


#include <math.h>

#include <mpi.h>

void parallel_precond_cg(Block *x, Block*b, MPI_Comm comm, MPI_Datatype vec);
void precond_conjugate_gradient(Grid *x, Grid *b);
void conjugate_gradient(Grid *x, Grid *b);

#endif
