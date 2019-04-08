#include "haloswaps.h"

void halo_swaps(Block* block, MPI_Comm comm, MPI_Datatype vecA, MPI_Datatype vecB, MPI_Datatype vecC, MPI_Datatype vecD, MPI_Datatype rowA, MPI_Datatype rowB, int oddr, int oddv)
{
	//MPI_Barrier(comm);
	//printf("at startline %d\n",block->id);
	//MPI_Barrier(comm);
	if(oddr == 0)
	{
		if(block->nbrup > -1 && block->nbrdown > -1)
		{
			// Send up, receive from below
			MPI_Sendrecv(block->data2d[0], 1, rowA, block->nbrup, 0, block->halos2d[3], 1, rowA, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);
			// Send down, receive from above
			MPI_Sendrecv(block->data2d[block->height-1], 1, rowA, block->nbrdown, 0, block->halos2d[0], 1, rowA, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrup > -1)
		{
			// Send UP
			MPI_Send(block->data2d[0], 1,rowA, block->nbrup, 0, comm);
			// Receive from above
			MPI_Recv(block->halos2d[0], 1, rowA, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrdown > -1)
		{
			// Receive from below
			MPI_Recv(block->halos2d[3], 1, rowA, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);
			// Send down
			MPI_Send(block->data2d[block->height-1], 1, rowA, block->nbrdown, 0, comm);
		}
	}
	else
	{
		if(block->nbrup > -1 && block->nbrdown > -1)
		{
			// Send up, receive from below
			MPI_Sendrecv(&block->data2d[0][1], 1, rowB, block->nbrup, 0, &block->halos2d[3][1], 1, rowB, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);
			// Send down, receive from above
			MPI_Sendrecv(&block->data2d[block->height-1][1], 1, rowB, block->nbrdown, 0, &block->halos2d[0][1], 1, rowB, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrup > -1)
		{
			// Send UP
			MPI_Send(&block->data2d[0][1], 1, rowB, block->nbrup, 0, comm);
			// Receive from above
			MPI_Recv(&block->halos2d[0][1], 1, rowB, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrdown > -1)
		{
			// Receive from below
			MPI_Recv(&block->halos2d[3][1], 1, rowB, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);
			// Send down
			MPI_Send(&block->data2d[block->height-1][1], 1, rowB, block->nbrdown, 0, comm);
		}
	}
	MPI_Barrier(comm);
	//printf("made halfway %d\n",block->id);
//fflush(stdout);
	MPI_Barrier(comm);

	if (oddv == 0)
	{
		if(block->nbrleft > -1 && block->nbrright > -1)
		{
			// Send left, receive from right
			MPI_Sendrecv(block->data2d[0], 1, vecA, block->nbrleft, 0, block->halos2d[2], 1, vecC, block->nbrright, 0, comm, MPI_STATUS_IGNORE);
			//printf("myid is %d i send reked %d %d\n",block->id,block->nbrleft, block->nbrright);
			// Send right, receive from left
			MPI_Sendrecv(&block->data2d[0][block->width-1], 1, vecA, block->nbrright, 0, block->halos2d[1], 1, vecC, block->nbrleft, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrleft > -1)
		{
			// Send left
			MPI_Send(block->data2d[0], 1, vecA, block->nbrleft, 0, comm);
			//printf("myid is %d i sENt to %d\n",block->id,block->nbrleft);
			// Receive from left
			MPI_Recv(block->halos2d[1], 1, vecC, block->nbrleft, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrright > -1)
		{
			// Receive from right
			MPI_Recv(block->halos2d[2], 1, vecC, block->nbrright, 0, comm, MPI_STATUS_IGNORE);
			//printf("myid is %d i red to %d NOW SEND TO same\n",block->id,block->nbrright);
			// Send right
			MPI_Send(&block->data2d[0][block->width-1], 1, vecC, block->nbrright, 0, comm);
			//printf("myid is %d i sent to %d \n",block->id,block->nbrright);
		}
	}
	else
	{
		//printf("This shouldn't print (yet)\n");
		if(block->nbrleft > -1 && block->nbrright > -1)
		{
			// Send left, receive from right
			MPI_Sendrecv(block->data2d[1], 1, vecB, block->nbrleft, 0, &block->halos2d[2][1], 1, vecD, block->nbrright, 0, comm, MPI_STATUS_IGNORE);
			// Send right, receive from left
			MPI_Sendrecv(&block->data2d[1][block->width-1], 1, vecB, block->nbrright, 0, &block->halos2d[1][1], 1, vecD, block->nbrleft, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrleft > -1)
		{
			// Send left
			MPI_Send(block->data2d[1], 1, vecB, block->nbrleft, 0, comm);
		//	printf("myid is %d i sent to %d\n",block->id,block->nbrleft);
			// Receive left
			MPI_Recv(&block->halos2d[1][1], 1, vecD, block->nbrleft, 0, comm, MPI_STATUS_IGNORE);
		}
		else if(block->nbrright > -1)
		{
			// Receive from right
			//printf("myid is %d i sent to %d\n",block->id,block->nbrright);
			MPI_Recv(&block->halos2d[2][1], 1, vecD, block->nbrright, 0, comm, MPI_STATUS_IGNORE);
			//printf("myid is %d\n",block->id);
			// Send right
			MPI_Send(&block->data2d[1][block->width-1], 1, vecD, block->nbrright, 0, comm);
		}
	}
}

void create_vec(Block *block, MPI_Datatype *vec)
{
	MPI_Type_vector(block->height, 1, block->width, MPI_DOUBLE, vec);
	MPI_Type_commit(vec);
}

void create_vecs(Block *block, MPI_Datatype *vecA, MPI_Datatype *vecB, MPI_Datatype *vecC, MPI_Datatype *vecD)
{
	if(block->height%2 == 0)
	{
  	MPI_Type_vector(block->height/2, 1, 2*block->width, MPI_DOUBLE, vecA);
  	MPI_Type_vector(block->height/2, 1, 2*block->width, MPI_DOUBLE, vecB);
  	MPI_Type_vector(block->height/2, 1, 2, MPI_DOUBLE, vecC);
  	MPI_Type_vector(block->height/2, 1, 2, MPI_DOUBLE, vecD);
		MPI_Type_commit(vecA);
		MPI_Type_commit(vecB);
		MPI_Type_commit(vecC);
		MPI_Type_commit(vecD);
	}
	else
	{
  	MPI_Type_vector((block->height/2)+1, 1, 2*block->width, MPI_DOUBLE, vecA);
  	MPI_Type_vector(block->height/2, 1, 2*block->width, MPI_DOUBLE, vecB);
  	MPI_Type_vector((block->height/2)+1, 1, 2, MPI_DOUBLE, vecC);
  	MPI_Type_vector(block->height/2, 1, 2, MPI_DOUBLE, vecD);
		MPI_Type_commit(vecA);
		MPI_Type_commit(vecB);
		MPI_Type_commit(vecC);
		MPI_Type_commit(vecD);
	}
}

void create_rows(Block *block, MPI_Datatype *rowA, MPI_Datatype *rowB)
{
	if(block->width%2 == 0)
	{
  	MPI_Type_vector(block->width/2, 1, 2, MPI_DOUBLE, rowA);
  	MPI_Type_vector(block->width/2, 1, 2, MPI_DOUBLE, rowB);
		MPI_Type_commit(rowA);
		MPI_Type_commit(rowB);
	}
	else
	{
  	MPI_Type_vector((block->height/2)+1, 1, 2, MPI_DOUBLE, rowA);
  	MPI_Type_vector(block->height/2, 1, 2, MPI_DOUBLE, rowB);
		MPI_Type_commit(rowA);
		MPI_Type_commit(rowB);
	}
}

void halo_swapping(Block* block, MPI_Comm comm, MPI_Datatype vec)
{	

	if(block->nbrup > -1 && block->nbrdown > -1)
	{
		// Send up, receive from below
		MPI_Sendrecv(block->data2d[0], block->width, MPI_DOUBLE, block->nbrup, 0, block->halos2d[3], block->width, MPI_DOUBLE, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);
		// Send down, receive from above
		MPI_Sendrecv(block->data2d[block->height-1], block->width, MPI_DOUBLE, block->nbrdown, 0, block->halos2d[0], block->width, MPI_DOUBLE, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
	}
	else if(block->nbrup > -1)
	{
		// Send UP
		MPI_Send(block->data2d[0], block->width, MPI_DOUBLE, block->nbrup, 0, comm);
		// Receive from above
		MPI_Recv(block->halos2d[0], block->width, MPI_DOUBLE, block->nbrup, 0, comm, MPI_STATUS_IGNORE);
	}
	else if(block->nbrdown > -1)
	{
		// Receive from below
		MPI_Recv(block->halos2d[3], block->width, MPI_DOUBLE, block->nbrdown, 0, comm, MPI_STATUS_IGNORE);
		// Send down
		MPI_Send(block->data2d[block->height-1], block->width, MPI_DOUBLE, block->nbrdown, 0, comm);
	}
	
	if(block->nbrleft > -1 && block->nbrright > -1)
	{
		// Send left, receive from right
		MPI_Sendrecv(block->data2d[0], 1, vec, block->nbrleft, 0, block->halos2d[2], block->height, MPI_DOUBLE, block->nbrright, 0, comm, MPI_STATUS_IGNORE);
		// Send right, receive from left
		MPI_Sendrecv(&block->data2d[0][block->width-1], 1, vec, block->nbrright, 0, block->halos2d[1], block->height, MPI_DOUBLE, block->nbrleft, 0, comm, MPI_STATUS_IGNORE);
	}
	else if(block->nbrleft > -1)
	{
		// Send left
		MPI_Send(block->data2d[0], 1, vec, block->nbrleft, 0, comm);
		// Receive from left
		MPI_Recv(block->halos2d[1], block->height, MPI_DOUBLE, block->nbrleft, 0, comm, MPI_STATUS_IGNORE);
	}
	else if(block->nbrright > -1)
	{
		// Receive from right
		MPI_Recv(block->halos2d[2], block->height, MPI_DOUBLE, block->nbrright, 0, comm, MPI_STATUS_IGNORE);
		//printf("myid is %d i red to %d NOW SEND TO same\n",block->id,block->nbrright);
		// Send right
		MPI_Send(&block->data2d[0][block->width-1], 1, vec, block->nbrright, 0, comm);
		//printf("myid is %d i sent to %d \n",block->id,block->nbrright);
	}
}
