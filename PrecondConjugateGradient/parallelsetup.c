#include "parallelsetup.h"

void alloc_block(Block* my_block, int size, int myid)
{
	my_block->id = myid;
	my_block->size = size;
	if(my_block->id == 2 )
	{
		my_block->width = size-1;
		my_block->height = size;
	}
	else if(my_block->id == 6 || my_block->id == 7)
	{
		my_block->width = size;
		my_block->height = size-1;
	}
	else
	{
		my_block->width = size;
		my_block->height = size;
	}
	find_my_neighbours(my_block);

	my_block->data1d = malloc(sizeof(double)*my_block->width*my_block->height);
	my_block->data2d = malloc(sizeof(double*)*(my_block->height));

	int halo_number[4] = {0,0,0,0};
	if( my_block->nbrup != -1 )
	{
		halo_number[0]++;
	}
	if (my_block->nbrdown != -1)
	{
		halo_number[3]++;
	}	
	if(	my_block->nbrright != -1)
	{
		halo_number[2]++;
	}
	if ( my_block->nbrleft != -1)
	{
		halo_number[1]++;
	}

	my_block->halos1d = malloc(sizeof(double)*(((halo_number[1]+halo_number[2])*my_block->height)+((halo_number[3]+halo_number[0])*my_block->width)));
	my_block->halos2d = malloc(sizeof(double*)*(4));

	int i;
	for ( i=0; i<my_block->height; i++)
	{
		my_block->data2d[i] = &my_block->data1d[i*my_block->width];
	}
	int counter=0;
	for (i=0; i<4; i++)
	{
		my_block->halos2d[i] = &my_block->halos1d[counter];
		if(halo_number[i]==1)
		{
			if(i==0)
			{
				counter+=my_block->width;
			}		
			else
			{
				counter+=my_block->height;
			}
		}	
	}
	my_block->halo_size = counter;
}

void initiate_block(Block* block, MPI_Comm comm)
{
	int i;
	MPI_Barrier(comm);
	for( i=0; i<block->height*block->width; i++)
	{
		block->data1d[i] = 0.0;
	}
	if(block->nbrup != -1)
	{
		for(i=0;i<block->width;i++)
		{
			block->halos2d[0][i] = 0.0;
		}
	}
	if(block->nbrdown != -1)
	{
		for(i=0;i<block->width;i++)
		{
			block->halos2d[3][i] = 0.0;
		}
	}
	if(block->nbrleft != -1)
	{
		for(i=0;i<block->height;i++)
		{
			block->halos2d[1][i] = 0.0;
		}
	}
	if(block->nbrright != -1)
	{
		block->halos2d[2][i] = 0.0;
	}
}

void fill_blockB(Block* block) 
{
	int i;
	if (block->id != 2)
	{
		for( i=0; i<block->height*block->width; i++)
		{
			block->data1d[i] = 0.0;
		}
	}
	else
	{
		for( i=0; i<block->height*block->width; i++)
		{
			block->data1d[i] = 0.0;
		}
		for(i=block->width-1; i<block->height*block->width; i+=block->width)
		{
			block->data1d[i] = 1.0;
		}
	}
}

void print_block(Block *block, MPI_Comm comm)
{
	MPI_Barrier(comm);
	int i,j;
	if (block->id == 0)
	{
		int wid,hei;
		Grid grid;
		alloc_serial_grid(&grid, block->size);
		for (i=0; i<block->height; i++)
		{
			for(j=0; j<block->width; j++)
			{
				grid.data2d[i][j] = block->data2d[i][j];
			}
		}
		for(i=1; i<8; i++)
		{
			MPI_Recv(&hei, 1, MPI_INT, i, 1, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&wid, 1, MPI_INT, i, 2, comm, MPI_STATUS_IGNORE);
			if(i<3)
			{
				for(j=0; j<hei; j++)
				{
					MPI_Recv(&grid.data2d[j][grid.size*i], wid, MPI_DOUBLE, i, 3, comm, MPI_STATUS_IGNORE);
				}
			}
			else if(i<4)
			{
				for(j=0; j<hei; j++)
				{
					MPI_Recv(grid.data2d[block->size+j], wid, MPI_DOUBLE, i, 3, comm, MPI_STATUS_IGNORE);
				}

			}
			else if(i<6)
			{
				for(j=0; j<hei; j++)
				{
					MPI_Recv(&grid.data2d[(2*block->size)+j][grid.size*(i-4)], wid, MPI_DOUBLE, i, 3, comm, MPI_STATUS_IGNORE);
				}
			}
			else
			{
				for(j=0; j<hei; j++)
				{
					MPI_Recv(&grid.data2d[(3*grid.size)+j][grid.size*(i-6)], wid, MPI_DOUBLE, i, 3, comm, MPI_STATUS_IGNORE);
				}
			}
		}
		print_full_grid(&grid);
		free_grid(&grid);
	}
	else
	{
		for(i=1; i<8; i++)
		{
			if(block->id == i)
			{
				MPI_Ssend(&block->height, 1, MPI_INT, 0, 1, comm);
				MPI_Send(&block->width, 1, MPI_INT, 0, 2, comm);
				for(j=0; j<block->height ; j++)
				{
					MPI_Send(block->data2d[j], block->width, MPI_DOUBLE, 0, 3, comm);
				}
			}
		}	
	}
}

void free_block(Block* my_block)
{
	free(my_block->data1d);
	free(my_block->data2d);
	free(my_block->halos1d);
	free(my_block->halos2d); 
}

void find_my_neighbours(Block* my_block)
{
	if(my_block->id == 0)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = 3;	
		my_block->nbrright = 1;
		my_block->nbrleft = -1;

	}
	if(my_block->id == 1)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = -1;	
		my_block->nbrright = 2;
		my_block->nbrleft = 0;

	}	
	if(my_block->id == 2)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = -1;	
		my_block->nbrright = -2;
		my_block->nbrleft = 1;

	}	
	if(my_block->id == 3)
	{
		my_block->nbrup = 0;
		my_block->nbrdown = 4;	
		my_block->nbrright = -1;
		my_block->nbrleft = -1;

	}	
	if(my_block->id == 4)
	{
		my_block->nbrup = 3;
		my_block->nbrdown = 6;	
		my_block->nbrright = 5;
		my_block->nbrleft = -1;

	}	
	if(my_block->id == 5)
	{
		my_block->nbrup = -1;
		my_block->nbrdown = 7;	
		my_block->nbrright = -1;
		my_block->nbrleft = 4;

	}	
	if(my_block->id == 6)
	{
		my_block->nbrup = 4;
		my_block->nbrdown = -2;	
		my_block->nbrright = 7;
		my_block->nbrleft = -1;

	}	
	if(my_block->id == 7)
	{
		my_block->nbrup = 5;
		my_block->nbrdown = -2;	
		my_block->nbrright = -1;
		my_block->nbrleft = 6;

	}	
}



double pa_nbr_up(Block *b, int x, int y)
{
	if(x == 0)
	{
		if(b->nbrup == -1)
		{
			return b->data2d[x+1][y];
		}
		else
		{
			return b->halos2d[0][y];
		}
	}
	return b->data2d[x-1][y];
}

double pa_nbr_down(Block *b, int x, int y)
{
	if (x==(b->height-1))
	{
		if(b->nbrdown == -1)
		{
			return b->data2d[x-1][y];
		}
		else
		{
			return b->halos2d[3][y];
		}
	}
	return b->data2d[x+1][y];
}

double pa_nbr_left(Block *b, int x, int y)
{
	if(y==0)
	{
		if(b->nbrleft == -1)
		{
			return b->data2d[x][y+1];
		}
		else
		{
			return b->halos2d[1][x];
		}
	}
	return b->data2d[x][y-1];
}

double pa_nbr_right(Block *b, int x, int y)
{
	if (y==b->width-1)
	{
		if(b->nbrright == -1)
		{
			return b->data2d[x][y-1];
		}
		else
		{
			return b->halos2d[2][x];
		}
	}
	return b->data2d[x][y+1];
}




