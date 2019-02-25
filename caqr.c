#include "caqr.h"
#include "hhorth.h"
#include "parallelsetup.h"
#include "tsqr.h"

#include <string.h>

#include <mpi.h>

void caqr(DenseMatrix* Rfin, DenseMatrix *Qfin, int nprocs, int myid, int b)
/* Function to apply communication avoiding qr factorization to a matrix stored
in Rfin across nprocs number of processors. The values of Rfin will be updated
with the values of R for the R factorization. */
{
	int s,e;
	int i,j,k,p;
	int os, oe;

	//Calculation of how far into communication step the processor must go:
	int comm_steps = int_log2(nprocs);
	int mylevel = my_steps(nprocs, myid, comm_steps);

	// Matrices declared for implicit calculation of Q and R.
	DenseMatrix Q[mylevel];
	DenseMatrix R[mylevel];


	for(j=0; j<Rfin->nColumns/b; j++)
	{
		//Each processor takes it's chunk of matrix.
		decomp1d(j*b, Rfin->nRows, nprocs, myid, &s, &e);

	// I created a deep copy matrix W so that each block would be contiguously
	// stored in memory. With more time I would test if this was worthwhile or if
	// it would be faster to deal directly with the Rfin matrix.
		DenseMatrix W = CreateNullMatrix(1+e-s,b);
		for (i=0; i<b; i++)
		{
			memcpy(W.entry[i],&Rfin->entry[(j*b)+i][s],sizeof(double)*(1+e-s));
		}

		// The first iterations are done outside of a loop.
		R[0] = CreateNullMatrix(1+e-s,b);
		Q[0] = CreateNullMatrix(1+e-s,b);
		mod_hhorth(&W,&Q[0],&R[0]);
		// We now iterate through, only saving the final R matrix
		// but saving all q representatives.
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

		// Now R has been calculated we update the trailing matrix
		for (p=j+1;p<Rfin->nColumns/b;p++)
		{
			for (i=0; i<b; i++)
			{
				memcpy(W.entry[i],&Rfin->entry[(p*b)+i][s],sizeof(double)*(1+e-s));
			}
			caQtA(Q, &W, myid, comm_steps, mylevel, b, j,s);
			for (i = 0; i<b; i++)
			{
				memcpy(&Rfin->entry[(p*b)+i][s],W.entry[i],sizeof(double)*W.nRows);
			}
		}

		// We store the final R matrix before freeing:
		if (myid == 0)
		{
			for (i=0; i<b; i++)
			{
				memcpy(Rfin->entry[(j*b)+i],R[mylevel-1].entry[i],sizeof(double)*(i+1));
			}
			FreeMatrix(&R[mylevel-1]);
		}

		// Needs to be changed, at the moment exchanges all updates of Rfin,
		// Should only do necessary ones.
		for (p=0;p<nprocs;p++)
		{
			os = s;
			oe = e;
			MPI_Bcast(&os, 1, MPI_INT, p,MPI_COMM_WORLD);
			MPI_Bcast(&oe, 1, MPI_INT, p,MPI_COMM_WORLD);
			for (k=0;k<Rfin->nColumns;k++)
			{
				MPI_Bcast(&Rfin->entry[k][os], 1+oe-os, MPI_DOUBLE, p,MPI_COMM_WORLD);
			}
		}


		// The last dynamically allocated memory freed before finish of iteration.
		for (i = 0; i < mylevel; i++)
		{
			FreeMatrix(&Q[i]);
		}
		FreeMatrix(&W);
	}

	// Can be uncommented for error checking:
	// if (myid == 0)
	// {
	// 	printf("Rfin[10][%d] = %lf\n",0,Rfin->entry[0][0]);
	// 	for (i=0;i< 15;i++){
	// 		for (j=0; j< 15; j++){
	// 			printf("  %lf  ",Rfin->entry[j][i]);
	// 		} printf("\n");}
	// }
}




void caQtA(DenseMatrix Q[], DenseMatrix *W, int myid, int comm_steps,int mylevel, int b, int j, int s)
/* Function for applying implicitly stored Q^T across trailing matrix for
communication avoiding qr factorization. */
{

	int i,k;
	DenseMatrix Bbar;

	// First application of Q^T across trailing matrix.
	Apply_QT(W,&Q[0]);

	// Bbar used for combining higher level Q matrices
	//if (mylevel>1){}
	Bbar = CreateNullMatrix(2*b,b);
	for ( i=0;i<b;i++)
	{
		memcpy(Bbar.entry[i],W->entry[i],sizeof(double)*b);
	}

	// Enter communication levels again for application of higher level Q matrices.
	for (k=1; k<=comm_steps; k++ )
	{
		if (k < mylevel)
		{
			for(i=0;i<b;i++)
			{
				MPI_Recv(&Bbar.entry[i][b], b, MPI_DOUBLE, myid+int_pow(2,k-1), 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			ApplyQT(&Bbar,&Q[k],&Bbar);
			for(i=0;i<b;i++)
			{
				memcpy(W->entry[i],Bbar.entry[i],sizeof(double)*b);
				MPI_Send(&Bbar.entry[i][b], b, MPI_DOUBLE, myid+int_pow(2,k-1), 8, MPI_COMM_WORLD);
			}
		}
		else if (k == mylevel)
		{
			for(i=0;i<b;i++)
			{
				MPI_Send(W->entry[i],b,MPI_DOUBLE,myid-int_pow(2,k-1),7,MPI_COMM_WORLD);
			}
			for(i=0;i<b;i++)
			{
				MPI_Recv(W->entry[i], b, MPI_DOUBLE, myid-int_pow(2,k-1), 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}

	FreeMatrix(&Bbar);
}

void Apply_QT(DenseMatrix* B,DenseMatrix* Q)
/* Function to apply a Q^T as calculated for an nxb block
to each other block of a matrix*/
{
	int i,j;
	for (i=0; i<B->nColumns;i++)
	{
		for (j=0; j<Q->nColumns; j++)
		{
			Pk_TimesX(B->entry[i],Q->entry[j],B->nRows,j);
		}
	}
}
