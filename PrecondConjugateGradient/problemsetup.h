#ifndef PROBLEMSETUP_H
#define PROBLEMSETUP_H

#include <stdlib.h>
#include <stdio.h>

#include "structs.h"

/*	problemsetup.h -- header file for problemsetup.c
 *
 *	Author: John Cormican
 *
 */


void alloc_serial_grid(Grid *grid, int size);
void fill_grid(Grid *grid);
void define_b(Grid* x, Grid* b);

double nbr_up(Grid *grid, int x, int y);
double nbr_down(Grid *grid, int x, int y);
double nbr_left(Grid *grid, int x, int y);
double nbr_right(Grid *grid, int x, int y);


void print_full_grid(Grid* grid);
void free_grid(Grid *grid);

#endif
