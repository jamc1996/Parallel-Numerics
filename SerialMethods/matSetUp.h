#ifndef MATSETUP_H
#define MATSETUP_H

#include "dense.h"

// Functions to create matrices
DenseMatrix CreateRandomMatrix(int n, int m);
DenseMatrix CreateNullMatrix(int n, int m);
DenseMatrix CreateLT(int n, int m);

// Functions to update values of a matrix
void SetToZeros(DenseMatrix *A);
void SetToRandom(DenseMatrix *A);

// Used in error checking for small matrices.
int PrintMatrix(DenseMatrix *A);

// Memory freeing.
void FreeMatrix(DenseMatrix *A);

#endif
