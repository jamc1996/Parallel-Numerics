#include "vectorops.h"

void par_calc_error_vector(Block* r, Block* x, Block* b)
{
	int i,j;
	for( i=0; i<r->height; i++)
	{
		for (j=0; j<r->width; j++)
		{
			r->data2d[i][j] = b->data2d[i][j] - par_A_times_xj(x, i, j);
		}
	}
}

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

void par_find_coords(Block*r, int index, int* x1, int* x2)
{
	*x1 = index/r->width;
	*x2 = index%r->width;
}


void linear_transform(double alpha, Grid* x, double beta, Grid*p)
{
	int i;
	for( i=0; i<x->total_size; i++)
	{
		x->data1d[i] = alpha*x->data1d[i] + beta*p->data1d[i];
	}
}

void par_linear_transform(double alpha, Block* x, double beta, Block *p)
{
	int i,j;
	for( i=0; i<x->height; i++)
	{
		for (j=0; j<x->width;j++)
		{
			x->data2d[i][j] = alpha*x->data2d[i][j] + beta*p->data2d[i][j];
		}
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

void copy_block(Block *p, Block *r)
{
	int i;
	for( i=0; i<p->height*p->width; i++)
	{
		p->data1d[i] = r->data1d[i];
	}
	for( i=0; i<p->halo_size; i++)
	{
		p->halos1d[i] = r->halos1d[i];
	}
}

double A_times_xj(Grid *grid, int x, int y)
{
	return (4*grid->data2d[x][y])-(nbr_up(grid,x,y) + nbr_down(grid,x,y) + nbr_left(grid,x,y) + nbr_right(grid,x,y)); 
}

double par_A_times_xj(Block *grid, int x, int y)
{
	return (4*grid->data2d[x][y])-(pa_nbr_up(grid,x,y) + pa_nbr_down(grid,x,y) + pa_nbr_left(grid,x,y) + pa_nbr_right(grid,x,y)); 
}


double noDiag_times_xj(Grid *grid, int x, int y)
{
	return -(nbr_up(grid,x,y) + nbr_down(grid,x,y) + nbr_left(grid,x,y) + nbr_right(grid,x,y)); 
}

double par_noDiag_times_xj(Block *grid, int x, int y)
{
	return -(pa_nbr_up(grid,x,y) + pa_nbr_down(grid,x,y) + pa_nbr_left(grid,x,y) + pa_nbr_right(grid,x,y)); 
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

void par_A_mul(Block* Ap, Block* p)
{
	int i,j;
	for( i=0; i<p->height; i++)
	{
		for(j=0; j<p->width; j++)
		{
			Ap->data2d[i][j] = par_A_times_xj(p, i, j);
		}
	}
}

double par_inner_prod(Block* x, Block* y)
{
	int i;
	double inner_prod = 0;
	for(i =0; i<x->width*x->height; i++)
	{
		inner_prod+=x->data1d[i]*y->data1d[i];
	}
	return inner_prod;
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

