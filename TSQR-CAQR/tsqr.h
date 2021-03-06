#ifndef TSQR_H
#define TSQR_H

#include "dense.h"

/*	tsqr.h -- header file for tsqr.c
 *
 *	Author: John Cormican
 *
 * */

int int_pow(int base, int exp);
int int_log2(int nprocs);
int my_steps(int nprocs, int myid, int comm_steps);

void tsqr(DenseMatrix* Rfin, int nprocs, int myid, int b);
void QtA(DenseMatrix Q[],DenseMatrix* B, DenseMatrix *W, int myid, int comm_steps, int mylevel, int b);

#endif
