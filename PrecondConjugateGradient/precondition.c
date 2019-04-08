#include "precondition.h"
#include "vectorops.h"

#include <stdio.h>
#include <stdlib.h>

void par_rb_gauss_seidel(Block* p, Block* r, MPI_Comm comm, MPI_Datatype vec)
{

	double sigma;
	int i, j;
	int a;
	for(a=0; a<1000; a++)
	{
		//Red nodes
		for (i=0;i<r->height;i++)
		{
			for(j=0;j<r->width;j++)
			{
				if((i+j)%2==0)
				{
					sigma = par_noDiag_times_xj(p, i, j);
					p->data2d[i][j] = 0.25*(r->data2d[i][j] -sigma);
				}
			}
		}	
		//Black nodes
		for (i=0;i<r->height;i++)
		{
			for(j=0;j<r->width;j++)
			{
				if((i+j)%2==1)
				{
					sigma = par_noDiag_times_xj(p, i, j);
					p->data2d[i][j] = 0.25*(r->data2d[i][j] -sigma);
				}
			}
		}	
		halo_swapping(p, comm, vec);
	}

}


void rb_gauss_seidel(Grid* p, Grid* r)
{
	double sigma=1.0;
	int i;
	int x1, x2;
	int a;
	for(a=0;a<1000;a++)
	{
		//Red nodes
		for (i=0;i<r->total_size;i++)
		{
			find_coords(p, i, &x1, &x2);
			if((x1+x2)%2==0)
			{
				sigma = noDiag_times_xj(p, x1, x2);
				p->data1d[i] = 0.25*(r->data1d[i]-sigma);
			}
		}
		//Black nodes
		for (i=0;i<r->total_size;i++)
		{
			find_coords(p, i, &x1, &x2);
			if((x1+x2)%2==1)
			{
				sigma = noDiag_times_xj(p, x1, x2);
				p->data1d[i] = 0.25*(r->data1d[i]-sigma);
			}
		}
		
	}
	
}
