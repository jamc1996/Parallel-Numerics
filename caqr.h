#ifndef CAQR_H
#define CAQR_H

#include "dense.h"

void caqr(DenseMatrix* Rfin, DenseMatrix *Qfin, int nprocs, int myid, int b);
void caQtA(DenseMatrix Q[], DenseMatrix *W, int myid, int comm_steps,int mylevel, int b, int j, int s);
void Apply_QT(DenseMatrix* B,DenseMatrix* Q);

#endif
