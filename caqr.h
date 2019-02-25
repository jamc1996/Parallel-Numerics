#ifndef CAQR_H
#define CAQR_H

#include "dense.h"

#include <mpi.h>

void caqr(DenseMatrix* Rfin, int nprocs, int myid, int b, MPI_Comm comm, int nbrup, int nbrdown);
void caQtA(DenseMatrix Q[], DenseMatrix *W, int myid, int comm_steps,int mylevel, int b, int j, int s, MPI_Comm comm);
void Apply_QT(DenseMatrix* B,DenseMatrix* Q);

#endif
