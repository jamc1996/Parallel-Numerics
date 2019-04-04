#ifndef PROBLEMSETUP_H
#define PROBLEMSETUP_H

#include <stdlib.h>
#include <stdio.h>

typedef struct
{
	int num_sections;
	int len[4];
	int start[4];
	int size;
	int total_size;
	double * data1d;
	double ** data2d;

} Grid;



void alloc_serial_grid(Grid *grid, int size);
void define_b(Grid* x, Grid* b);
double nbr_up(Grid *grid, int x, int y);
double nbr_down(Grid *grid, int x, int y);
double nbr_left(Grid *grid, int x, int y);
double nbr_right(Grid *grid, int x, int y);


void fill_grid(Grid *grid);
void print_full_grid(Grid* grid);
void free_grid(Grid *grid);

#endif
