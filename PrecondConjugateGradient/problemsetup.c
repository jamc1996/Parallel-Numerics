#include "problemsetup.h"

void alloc_serial_grid(Grid *grid, int size)
{
	// Length of each row is number of grid points + ghost processors to left and right

	int i;	
	grid->num_sections = 4;
	grid->len[0] = (size*3)-1;
	grid->len[1] = size;
	grid->len[2] = size*2;
	grid->len[3] = size*2;
 

	grid->size = size;

	grid->total_size = size*(grid->len[0]+grid->len[1]+grid->len[2])+(size-1)*(grid->len[3]);

	grid->start[0] = 0;
	grid->start[1] = grid->size*grid->len[0];
	grid->start[2] = grid->start[1] + (grid->size*grid->len[1]);
	grid->start[3] = grid->start[2] + (grid->size*grid->len[2]);
	
	grid->data1d = malloc(sizeof(double)*grid->total_size);
	grid->data2d = malloc(sizeof(double*)*((size*4)-1));
	for ( i=0; i<size; i++)
	{
		grid->data2d[i] = &grid->data1d[i*grid->len[0]];
	}
	for (i=size; i<2*size; i++ )
	{
		grid->data2d[i] = &grid->data1d[grid->start[1]+((i-size)*grid->len[1])];
	}
	for (i=size*2;i<size*3;i++)
	{
		grid->data2d[i] = &grid->data1d[grid->start[2]+((i-size*2)*grid->len[2])];
	}
	for (i=size*3;i<size*4-1;i++)
	{
		grid->data2d[i] = &grid->data1d[grid->start[3]+((i-size*3)*grid->len[2])];
	}
}

void define_b(Grid* x, Grid* b)
/* Define the b vector*/
{
	int i,j,k,counter=0;
	//A->nRows = x->total_size;
	//A->rowInd = malloc(sizeof(int)*A->nRows);
	//A->nEntries = 0;
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
				//A->rowInd[counter] = A->nEntries;
				b->data1d[counter] = 0.0;
				//A->nEntries += 5;
				if(k==0)
				{
					//b->data1d[counter]+=x->data2d[(x->size*i)+j][k+1];
					//A->nEntries--;
				}
				if(k==x->len[i]-1)
				{
					if(i==0)
					{
						b->data1d[counter]+=1.0;
						//A->nEntries--;
					}
					else
					{
						//b->data1d[counter]+=x->data2d[(x->size*i)+j][k-1];
						//A->nEntries--;
					}
				}
				if(j==0)
				{
					if(i==0)
					{
						//b->data1d[counter]+= x->data2d[(i*x->size)+j+1][k];
						//A->nEntries--;
					}
					else if(k >= x->len[i-1])
					{
						//b->data1d[counter]+= x->data2d[(i*x->size)+j+1][k];
						//A->nEntries--;
					} 
				}
				if(j==x->size-1)
				{
					if(k >= x->len[i+1])
					{
						//b->data1d[counter]+= x->data2d[(i*x->size)+j-1][k];
						//A->nEntries--;
					} 
				}
				counter++;
			}
		}
	}
}

double nbr_up(Grid *grid, int x, int y)
{
	if(x%grid->size == 0 && x/grid->size == 0)
	{
		return grid->data2d[x+1][y];
	}
	else if(x%grid->size == 0)
	{
		if (grid->len[(x/grid->size)-1] <= y)
		{
			return grid->data2d[x+1][y];
		}
	}
	return grid->data2d[x-1][y];
}

double nbr_down(Grid *grid, int x, int y)
{
	int i = x/grid->size;
	if (x==(4*grid->size)-2)
	{
		return 0.0;
	}
	if(x%grid->size == grid->size-1)
	{
		if (grid->len[i+1] <= y)
		{
			return grid->data2d[x-1][y];
		}
	}
	return grid->data2d[x+1][y];
}

double nbr_left(Grid *grid, int x, int y)
{
	if(y==0)
	{
		return grid->data2d[x][y+1];
	}
	return grid->data2d[x][y-1];
}

double nbr_right(Grid *grid, int x, int y)
{
	int i = x/grid->size;
	if (i==0)
	{
		if(y==grid->len[i]-1)
		{
			return 0.0;
		}
	}
	else if(y==grid->len[i]-1)
	{
		return grid->data2d[x][y-1];
	}
	return grid->data2d[x][y+1];
}

void fill_grid(Grid *grid)
{
	int i;

	for (i=0; i<grid->total_size; i++)
	{
		grid->data1d[i] = 0.0;
	}
}

void print_full_grid(Grid* grid)
{
	int i,j,k;
	for (i=0; i<4; i++)
	{
		for (j=0; j<grid->size; j++)
		{
			if(i==3 && j==grid->size-1)
			{
				break;
			}
			for (k=0; k<grid->len[i]; k++)
			{
				if(grid->data2d[(i*grid->size)+j][k]<0)
				{
					printf(" %.5lf ",grid->data2d[(i*grid->size)+j][k]);
				}
				else
				{
					printf("  %.5lf ",grid->data2d[(i*grid->size)+j][k]);
				}			
			}
			if(i==0)
			{
				printf("  %.5lf ",1.0);
			}
			printf("\n");
		}
	}
	for (i=0; i<grid->size*2; i++)
	{
		printf("  %.5lf ",0.0);
	}
	printf("\n");
}


void free_grid(Grid *grid)
{	
	free(grid->data2d);
	free(grid->data1d);
}
