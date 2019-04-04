#include "vectorops.h"



void calc_error_vector(Grid* r, Grid* x, Grid* b)
{
	int i;
	int x1, x2;
	for( i=0; i<r->total_size; i++)
	{
		find_coords(r, i, &x1, &x2);
		r->data1d[i] = b->data1d[i] - A_times_xj(x, x1, x2);
	}
}

void find_coords(Grid*r, int index, int* x1, int* x2)
{
	int i;
	for (i=0;i<4;i++)
	{
		if(index < r->len[i]*r->size)
		{	
			*x2 = index%r->len[i];
			*x1 = r->size*i + index/r->len[i];
			break;
		}
		else
		{
			index-=(r->len[i]*r->size);
		}
	}
}

void linear_transform(double alpha, Grid* x, double beta, Grid*p)
{
	int i;
	for( i=0; i<x->total_size; i++)
	{
		x->data1d[i] = alpha*x->data1d[i] + beta*p->data1d[i];
	}
}


void copy_grid(Grid *p, Grid *r)
{
	int i;
	for( i=0; i<r->total_size; i++)
	{
		p->data1d[i] = r->data1d[i];
	}
}

double A_times_xj(Grid *grid, int x, int y)
{
	return (4*grid->data2d[x][y])-(nbr_up(grid,x,y) + nbr_down(grid,x,y) + nbr_left(grid,x,y) + nbr_right(grid,x,y)); 
}

double noDiag_times_xj(Grid *grid, int x, int y)
{
	return -(nbr_up(grid,x,y) + nbr_down(grid,x,y) + nbr_left(grid,x,y) + nbr_right(grid,x,y)); 
}

void A_mul(Grid* Ap, Grid* p)
{
	int i;
	int x1, x2;
	for( i=0; i<p->total_size; i++)
	{
		find_coords(p, i, &x1, &x2);
		Ap->data1d[i] = A_times_xj(p, x1, x2);
	}
}

double inner_prod(Grid* x, Grid* y)
{
	int i;
	double inner_prod = 0;
	for( i=0; i<x->total_size; i++)
	{
		inner_prod += y->data1d[i]*x->data1d[i];
	}
	return inner_prod;
}

