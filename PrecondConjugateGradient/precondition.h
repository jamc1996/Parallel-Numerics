#ifndef PRECONDITION_H
#define PRECONDITION_H
#include "problemsetup.h"
typedef struct
{
	int nRows;
	int nEntries;
	int *rowInd;
	int *colInd;
	double *entries;

} SparseMatrix;

double lower_A_times_xj(Grid* grid, int x, int y);
void allocate_L(SparseMatrix *L, Grid *x);
void fill_L(SparseMatrix *L, Grid *x);
void rb_gauss_seidel(Grid* p, Grid* r);
#endif
