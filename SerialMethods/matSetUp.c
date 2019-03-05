#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "matSetUp.h"
#include "dense.h"

/*
		matSetUp.c -- 	program provides functions for creating, reseting and freeing dense matrices.

		Author:			John Cormican

		Purpose:		To provide functions necessary for managing the DenseMatrix structures.

		Usage:			Methods called from all other programs depending on what matrices they require.
*/


DenseMatrix CreateRandomMatrix(int n, int m)
/* Function to allocate memory for and return a DenseMatrix filled with
random values between -20 and 20. */
{
  DenseMatrix A;
  int i;
  A.nRows = n;
  A.nColumns = m;
  A.data_ = malloc(sizeof(double)*n*m);//all entries stored in t.data_.
  for (i=0;i<n*m;i++)
  {
    A.data_[i]= (drand48()*40) - 20;
  }
  A.entry = malloc(sizeof(double*)*m);//space allocated for m columns.
  for (i=0;i<m;i++)
  {
    A.entry[i] = &A.data_[i*n];
  }
  return A;
}

DenseMatrix CreateNullMatrix(int n, int m)
/* Function to create an n x m matrix and fill it with zeros. */
{
  DenseMatrix A;
  int j;
  A.nRows = n;
  A.nColumns = m;
  //A.data_ is used to ensure contiguous block memory used.
  A.data_ = malloc(sizeof(double)*n*m); //all entries stored in t.data_.
  memset(A.data_,0.0,sizeof(double)*n*m);

  // A stored column major:
  A.entry = malloc(sizeof(double*)*m);//space allocated for m columns.
  for (j=0;j<m;j++)
  {
    A.entry[j] = &A.data_[j*n];
  }

  /* -------> Warning memory allocated should be free later.  <----- */
  return A;
}

DenseMatrix CreateLT(int n, int m)
/* Creates enough room to efficiently store non zero values of m vectors with
zeros above the diagonal. */
{
  // Storage space calculated
	if(n<m){
		printf("Major problem\n");
	}
  DenseMatrix A;
	int sum = m*(n-m);
  int i;
  for ( i = 0; i < m; i++) {
    sum += i;
  }
  A.entry = malloc(sizeof(double*)*n);
  A.data_ = malloc(sizeof(double)*sum);

  // Organisation of pointers calculated.
  int offset = 0;
  for ( i = 0; i < m; i++) {
    offset+=i;
    A.entry[i] = &A.data_[(i*n)-offset];
  }
  return A;
}

int PrintMatrix(DenseMatrix *A)
/* Function used in error checking to print small matrices.
Warning --> Will be slow for larger matrices. */
{
  // Because Matrix is column major, but print requires to row by row this
  // function will not be very efficient.

  if (A->nRows*A->nColumns > 400) {
    printf("gramS: matrix size too large to print\n" );
    return 1;
  }

  int i,j;
  printf("\n" );
  for ( i = 0; i < A->nRows; i++) {
    for ( j = 0; j < A->nColumns; j++) {
      if(A->entry[j][i] < 0){
        printf(" %01.4lf ",A->entry[j][i] );
      }else{
        printf(" %02.4lf ",A->entry[j][i] );
      }
    }
    printf("\n" );
  }
  printf("\n" );

  return 0;
}

void FreeMatrix(DenseMatrix *A)
/* The memory dynamically allocated in DenseMatrix is freed. */
{
  free(A->data_);
  free(A->entry);
}


void SetToZeros(DenseMatrix *A)
/* memset used to set all values of DenseMatrix A to zero. */
{
  memset(A->data_,0.0,sizeof(double)*A->nColumns*A->nRows);
}

void SetToRandom(DenseMatrix *A)
/* All elements of a DenseMatrix A set to random values between -20 and 20. */
{
  int i;
  for (i=0;i<A->nRows*A->nColumns ;i++)
  {
    A->data_[i]= (drand48()*40) - 20;
  }
}
