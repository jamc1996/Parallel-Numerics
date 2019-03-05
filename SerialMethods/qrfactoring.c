
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "matops.h"
#include "matSetUp.h"
#include "qrfactoring.h"
#include "dense.h"

/*
		qrfactoring.c -- 	program provides functions for serial QR factorisation by
											both basic Gram-Schmidt and Householder.

		Author:			John Cormican

		Purpouse:		To demonstrate qr factorization methods.

		Usage:			Methods gs and hhorth called from main.c to perform qr factorisation.
								hhorth calls the other functions below as required.
*/



void gs(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R)
/* Implements basic version of Gram-Schmidt algorithm on matrix A. Updates
values in dense matrices Q and R to give the QR factorization. */
{
  SetToZeros(R);

  //First column of R updated.
  int i, j, k;
  R->entry[0][0] = two_norm(A->entry[0],A->nRows);
  if (R->entry[0][0] == 0) {
    return;
  }

  // First column of Q updated.
  scalar_div(A->entry[0],R->entry[0][0],Q->entry[0],Q->nRows);

  // Other columns filled in.
  for ( j = 1; j < A->nColumns; j++)
  {
    // R column j entries updated to adjust for previous steps.
    for (i = 0; i < j; i++)
    {
      R->entry[j][i] = inner_prod(A->entry[j],Q->entry[i],R->nRows);
    }

    // Q column j updated accordingly
    memcpy(Q->entry[j], A->entry[j], A->nRows * sizeof(double));
    for ( i = 0; i < j; i++) {
        for (k = 0; k < Q->nRows; k++) {
          Q->entry[j][k] -= R->entry[j][i]*Q->entry[i][k];
        }
    }

    //R[j][j] updated and Q entries normalised.
    R->entry[j][j] = two_norm(Q->entry[j],Q->nRows);
    if (R->entry[j][j] == 0) {
      printf("Error zero\n" );
      return;
    }
    scalar_div(Q->entry[j],R->entry[j][j],Q->entry[j],Q->nRows);
  }
}



void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R)
/* Implements QR factorization of matrix A via Householder orthogonalisation.
Updates the values of dense matrices Q and R. */
{
  // Iterators declared and place to store omega vectors. All omega vectors were
  // saved together as a lower triangular Matrix to save space.
  int i, k;
  DenseMatrix Omega = CreateLT(A->nRows, A->nColumns);

  for ( k = 0; k < A->nColumns; k++)
  {
    // Values of A_k are copied into R_k
    memcpy(R->entry[k],A->entry[k],(A->nRows)*sizeof(double));

    // Loop takes into account the effects of previous P_k transformations.
    if (k>0)
    {
        for (i = 0; i < k; i++)
        {
          Pk_TimesX(R->entry[k],Omega.entry[i],R->nRows,i);
        }
    }

    // Omega calulated to zero subdiagonal entries of k-th column of R.
    Calculate_Omega(Omega.entry[k],R->entry[k],k,R->nRows);

    // R_k updated according to P_k R_k
    Pk_TimesX(R->entry[k],Omega.entry[k],R->nRows,k);

    //Q columns calculated
    Q->entry[k][k]++;
    for ( i = k; i > -1; i--)
    {
      Pk_TimesX(Q->entry[k],Omega.entry[i],Q->nRows,i);
    }
  }

  FreeMatrix(&Omega);
}

void Calculate_Omega(Vector omeg_k, Vector A_k, int k, int n)
/* Function to calculate the omega_k vector needed
for P_k matrix for householder algorithm. */
{
  int i;

  // Values of omega updated
  double beta = two_norm(&A_k[k],n-k);
  if (A_k[k] < 0)
  {
    beta *= -1;
  }
  omeg_k[k] = beta+A_k[k];
  for (i = k+1; i < n; i++)
  {
    omeg_k[i]= A_k[i];
  }

  //omega normalised
  double norm = two_norm(&omeg_k[k], n-k);
  for ( i = k; i < n; i++) {
    omeg_k[i] /= norm;
  }
}

void Pk_TimesX(Vector rk, Vector omeg_k, int n, int k)
/* Function to perform multiplication of a vector by matrix P_k using
the omega_k vectors. */
{
  int i;

  // v = 2 r_k ^ T omega calculated
  double v = 0;
  for ( i = k; i < n; i++)
  {
    v += 2*rk[i]*omeg_k[i];
  }

  // r_k - omega_k v calculated.
  for ( i = k; i < n; i++)
  {
    rk[i] -= v*omeg_k[i];
  }
}
