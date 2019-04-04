#ifndef VECTOROPS_H
#define VECTOROPS_H

#include <stdio.h>

#include "problemsetup.h"

void calc_error_vector(Grid* r, Grid* x, Grid *b);
void find_coords(Grid*r, int index, int* x1, int* x2);

double A_times_xj(Grid *grid, int x, int y);

double inner_prod(Grid* x, Grid* y);

double A_times_xj(Grid *grid, int x, int y);
double precond_A_times_xj(Grid *grid, int x, int y);

void copy_grid(Grid *p, Grid *r);
double noDiag_times_xj(Grid *grid, int x, int y);;
void linear_transform(double alpha, Grid* x, double beta, Grid*p);
void A_mul(Grid* Ap, Grid* p);
#endif
