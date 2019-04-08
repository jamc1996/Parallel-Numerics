#ifndef VECTOROPS_H
#define VECTOROPS_H

#include <stdio.h>

#include "structs.h"
#include "parallelsetup.h"
#include "problemsetup.h"

void calc_error_vector(Grid* r, Grid* x, Grid *b);
void par_calc_error_vector(Block* r, Block* x, Block* b);

void find_coords(Grid*r, int index, int* x1, int* x2);
void par_find_coords(Block*r, int index, int* x1, int* x2);

double A_times_xj(Grid *grid, int x, int y);
double par_A_times_xj(Block *block, int x, int y);

double inner_prod(Grid* x, Grid* y);
double par_inner_prod(Block* x, Block* y);


void copy_grid(Grid *p, Grid *r);
void copy_block(Block *p, Block *r);

double noDiag_times_xj(Grid *grid, int x, int y);
double par_noDiag_times_xj(Block *grid, int x, int y);

void par_linear_transform(double alpha, Block* x, double beta, Block *p);
void linear_transform(double alpha, Grid* x, double beta, Grid*p);

void par_A_mul(Block* Ap, Block* p);
void A_mul(Grid* Ap, Grid* p);

#endif
