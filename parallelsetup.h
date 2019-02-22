#ifndef PARALLELSETUP_H
#define PARALLELSETUP_H

#include "dense.h"

int decomp1d(int offset, int n, int nprocs, int myid, int *sp, int *ep);
//void QtA(DenseMatrix* Q, DenseMatrix* A, int s, int step, int b);
#endif
