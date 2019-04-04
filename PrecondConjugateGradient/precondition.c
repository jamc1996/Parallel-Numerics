#include "precondition.h"
#include "vectorops.h"

#include <stdio.h>
#include <stdlib.h>

double lower_A_times_xj(Grid* grid, int x, int y)
{
	if(x>grid->size)
	{
		if(y==grid->len[x/grid->size] - 1)
		{
			return (4*grid->data2d[x][y]) + nbr_down(grid,x,y) + 2*nbr_left(grid,x,y); 
		}
	}

	return (4*grid->data2d[x][y]) + nbr_down(grid,x,y) + nbr_left(grid,x,y); 
}


void rb_gauss_seidel(Grid* p, Grid* r)
{
	double sigma=1.0;
	int i, j;
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
				sigma=0.0;
				
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
				sigma=0.0;

				sigma = noDiag_times_xj(p, x1, x2);
			
				p->data1d[i] = 0.25*(r->data1d[i]-sigma);
			}
		}
		
	}
	
}



void gauss_seidel(Grid* grid, Grid* b)
{
	double sigma;
	int i, j;
	int x1, x2;
	while(1==0)
	{
		for (i=0;i<grid->total_size;i++)
		{
			sigma=0;
			for (j=0; j<grid->total_size; j++)
			{
				find_coords(grid, j, &x1, &x2);
				sigma += A_times_xj(grid, x1, x2);
			}
			grid->data1d[i] = 0.25*(b->data1d[i]-sigma);
		}
		
	}
	
}



void allocate_L(SparseMatrix *L, Grid *x)
{
	L->nRows=x->total_size;
	L->rowInd = (int*)malloc(sizeof(int)*L->nRows);
	L->nEntries = (x->total_size*3) - (((x->size*4)-1) + ((x->size*3)-1) + x->size*1);
	L->entries = malloc(sizeof(double)*L->nEntries);
	L->colInd = (int*)malloc(sizeof(int)*L->nEntries);	

	int i;
	for(i=0; i<L->nEntries; i++)
	{
		printf("%lf %d\n",L->entries[i],L->colInd[i]);
	}
	fill_L(L,x);
}

void fill_L(SparseMatrix* L, Grid* x)
{
	int i,j, k;
	int counter=0;
	int entry_counter=0;
	for (i=0; i<4; i++)
	{
		for (j=0; j<x->size; j++)
		{
			if(i==3 && j==x->size-1)
			{
				break;
			}
			for (k=0; k<x->len[i]; k++)
			{
				//printf("%d, %d\n",counter, entry_counter);
				L->rowInd[counter] = entry_counter;
				if(k==0)
				{
					if (counter == 0)
					{
				//printf("is this breaking?\n"); 
						L->entries[entry_counter] = 4.0;
				//printf("is this breaking?\n"); 
						L->colInd[entry_counter] = counter; 
						entry_counter++;
						counter++;
						continue;
					}
					else
					{
						L->entries[entry_counter] = -1.0;
						if(j==0)
						{
							L->colInd[entry_counter] = counter-x->len[i-1];
						}
						else
						{
							L->colInd[entry_counter] = counter-x->len[i];
						}
						L->entries[entry_counter+1] = 4.0;
						L->colInd[entry_counter+1] = counter;
						entry_counter+=2;
						counter++;
						continue;
					}
				}
	
				if((i==0 && j==0) || (i==2 && j==0 && k>=x->len[i-1]) )
				{
						if(i==2 && k==x->len[i]-1)
						{
							L->entries[entry_counter] = -2.0;
							L->colInd[entry_counter] = counter-1;
							L->entries[entry_counter+1] = 4.0;
							L->colInd[entry_counter+1] = counter;
							entry_counter+=2;
							counter++;
							continue;
						}
						else
						{
							L->entries[entry_counter] = -1.0;
							L->colInd[entry_counter] = counter-1;
							L->entries[entry_counter+1] = 4.0;
							L->colInd[entry_counter+1] = counter;
							entry_counter+=2;
							counter++;
							continue;
						}
				}
					
				if(i==0 && j==x->size-1 && k>= x->len[1])
				{
						L->entries[entry_counter] = -2.0;
						L->colInd[entry_counter] = counter-x->len[i];
						L->entries[entry_counter+1] = -1.0;
						L->colInd[entry_counter+1] = counter -1;
						L->entries[entry_counter+2] = 4.0;
						L->colInd[entry_counter+2] = counter;
						entry_counter+=3;
						counter++;
						continue;
				}

				if(i>0 && j==x->len[i]-1)
				{
					L->entries[entry_counter] = -1;
					if(j==0)
					{
						L->colInd[entry_counter] = counter-x->len[i-1];
					}
					else
					{
						L->colInd[entry_counter] = counter-x->len[i];
					}
					L->entries[entry_counter+1] = -2.0;
					L->colInd[entry_counter+1] = counter-1;
					L->entries[entry_counter+2] = 4.0;
					L->colInd[entry_counter+2] = counter;
					entry_counter+=3;
					counter++;
					continue;
				}
				if(counter < x->len[0])
				{
					L->entries[entry_counter] = -1.0;
					L->colInd[entry_counter] = counter-1; 
					L->entries[entry_counter] = 4.0;
					L->colInd[entry_counter+1] = counter;
					entry_counter+=2;
					counter++;
					continue;
				}

				L->entries[entry_counter] = -1.0;
				if(j==0)
				{
					L->colInd[entry_counter] = counter-x->len[i-1];
				}
				else
				{
					L->colInd[entry_counter] = counter-x->len[i];
				}
				L->entries[entry_counter+1] = -1.0;
				L->colInd[entry_counter+1] = counter -1;
				L->entries[entry_counter+2] = 4.0;
				L->colInd[entry_counter+2] = counter;
				entry_counter+=3;
				counter++;
			}
		}
	}
}
