#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#include "hhorth.h"
#include "parallelsetup.h"
#include "dense.h"
#include "tsqr.h"



void QtA(DenseMatrix Q[],DenseMatrix* B, DenseMatrix *W, int myid, int comm_steps, int mylevel, int b)
{
	
	int i,k;
	ApplyQT(B,&Q[0],W);
	
	
	DenseMatrix Bbar = CreateNullMatrix(2*b,b);
	for ( i=0;i<b;i++)
	{
		memcpy(Bbar.entry[i],B->entry[i],sizeof(double)*b);
	}
	
	for (k=1; k<=comm_steps; k++ )
	{
		if (k < mylevel)
		{
			for(i=0;i<b;i++)
			{
				//printf("My id is %d, k is %d and i is %d\n",myid,k,i);
				MPI_Recv(&Bbar.entry[i][b], b, MPI_DOUBLE, myid+int_pow(2,k-1), 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			ApplyQT(&Bbar,&Q[k],&Bbar);			
			for(i=0;i<b;i++)
			{
				//printf("My id is %d, k is %d and i is %d\n",myid,k,i);
				MPI_Send(&Bbar.entry[i][b], b, MPI_DOUBLE, myid+int_pow(2,k-1), 8, MPI_COMM_WORLD);
				memcpy(B->entry[i],Bbar.entry[i],sizeof(double)*b);
			}
		}
		else if (k == mylevel)
		{
			for(i=0;i<b;i++)
			{
				//printf("My id is %d, k is %d and i is %d\n",myid,k,i);
				MPI_Send(B->entry[i],b,MPI_DOUBLE,myid-int_pow(2,k-1),7,MPI_COMM_WORLD);
			}
			for(i=0;i<b;i++)
			{
				//printf("My id is %d, k is %d and i is %d\n",myid,k,i);
				MPI_Recv(B->entry[i], b, MPI_DOUBLE, myid-int_pow(2,k-1), 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

		}
		//break;
	}


}

void tsqr(DenseMatrix* Rfin, DenseMatrix *Qfin, int nprocs, int myid, int b)
{
	int s,e;
	int i,k;

	//Each processor takes it's chunk of matrix.
	decomp1d(0, Rfin->nRows, nprocs, myid, &s, &e);
	DenseMatrix W = CreateNullMatrix(1+e-s,b);
	for (i=0; i<b; i++)
	{
		memcpy(W.entry[i],&Rfin->entry[i][s],sizeof(double)*(1+e-s));
	}


	int comm_steps = int_log2(nprocs);
	int mylevel = my_steps(nprocs, myid, comm_steps);
	DenseMatrix Q[mylevel];
	DenseMatrix R[mylevel];

	// First iterations are done
	R[0] = CreateNullMatrix(1+e-s,b);
	Q[0] = CreateNullMatrix(1+e-s,b);
	mod_hhorth(&W,&Q[0],&R[0]);

	// We now iterate through, only saving the final R matrix but saving all q representatives.
	for (i = 1; i<= comm_steps; i++)
	{
		if (i < mylevel)
			{
				for(k=0;k<b;k++)
				{
					MPI_Recv(&R[i-1].entry[k][b], k+1, MPI_DOUBLE, myid+int_pow(2,i-1), myid+int_pow(2,i-1), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				Q[i] = CreateNullMatrix(2*b,b);
				R[i] = CreateNullMatrix(2*b,b);
				mod_hhorth(&R[i-1],&Q[i],&R[i]);
				FreeMatrix(&R[i-1]);			
			}
			else if (i == mylevel)
			{
				for(k=0;k<b;k++)
				{
					MPI_Send(R[i-1].entry[k],k+1,MPI_DOUBLE,myid-int_pow(2,i-1),myid,MPI_COMM_WORLD);
				}
				FreeMatrix(&R[i-1]);
			}
	}
	
	DenseMatrix B;
	B = CreateNullMatrix(1+e-s,b);
	QtA(Q, &B, &W, myid, comm_steps, mylevel, b);
	if (myid == 0)
	{
		for ( i=0; i<b; i++ )
		{
			//memcpy(Rfin->entry[i],R[mylevel].entry[i],sizeof(double)*R[mylevel].nRows);
		}
	}
	else
	{
		for (i=0; i<b; i++)
		{
			//memset(&Rfin->entry[i][s],0.0,sizeof(double)*(1+e-s));
		}
	}
//	if (myid ==0){PrintMatrix(Rfin);}
}


int int_pow(int base, int exp)
{
	int i, res = 1;
	for (i=0;i<exp;i++)
	{
		res*=base;
	}
	return res;
}

int my_steps(int nprocs, int myid, int comm_steps)
{
	int i, count = 1;
	for (i=1;i<comm_steps+1;i++)
	{
		if (myid%int_pow(2,i)!=0)
		{
			return count;
		}
		count++;
	}
	return count;
}

int int_log2(int nprocs)
{
	int count = 0;
	while (nprocs>1)
	{
		nprocs/=2;
		count++;		
	}
	return count;
}
