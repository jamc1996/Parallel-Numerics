#include <math.h>
#include <stdio.h>

#include "matSetUp.h"
#include "matops.h"
#include "dense.h"

/*
		matSetUp.c -- 	program provides functions for matrix/vector operations.

		Author:			John Cormican

		Purpose:		To provide functions necessary the operations of QR factorisation.

		Usage:			Methods mostly called from qrfactoring.c to perform required
                matrix operations.
*/

double inner_prod(Vector Aj, Vector Qi, int n)
/* Calculates inner product of vector Aj and Qi of length n. */
{
  int k;
  double sum = 0;
  for ( k = 0; k < n; k++) {
    sum+= Aj[k]*Qi[k];
  }
  return sum;
}

double two_norm(Vector Aj,int m)
/* Calculates 2-norm of vector Aj of length n. */
{
  int i;
  double sum = 0;
  for ( i = 0; i < m; i++) {
    sum+= Aj[i]*Aj[i];
  }
  return sqrt(sum);
}

void scalar_div(Vector Ai,double Rii, Vector Qj,int j)
/* Function to divide all entries in a vector by a scalar */
{
  int i;
  for ( i = 0; i < j; i++) {
    Qj[i] = Ai[i]/Rii;
  }
}

void scalar_mult(Vector Ai,double Rii, Vector Qj,int j)
/* Function to multiply all entries in a vector by a scalar */
{
  int i;
  for ( i = 0; i < j; i++) {
    Qj[i] = Ai[i]*Rii;
  }
}

void MatrixMultiply(DenseMatrix *Q, DenseMatrix *R, DenseMatrix *QR)
/* Function multiplies matrices Q and R and strores result in QR. */
{
  int i,j,k;

  // This function is also quite slow, but is only used in error checking.
  // I would have transposed Q to improve cache performance but this would
  // have involved a lot of restructuring.
  for ( j = 0; j < QR->nColumns; j++) {
    for ( i = 0; i < QR->nRows ; i++) {
      for ( k = 0; k < Q->nColumns; k++) {
          QR->entry[j][i] += Q->entry[k][i]*R->entry[j][k];
      }
    }
  }
  return;
}

double MaxMatValue(DenseMatrix *A)
/* Function to find the maximum absolute value of an entry in DenseMatrix A. */
{
	int i,j;
	double max = 0;

  // Values are iterated through columnwise
	for (j = 0; j<A->nColumns; j++){
		for (i = 0; i< A->nRows; i++){
			if (A->entry[j][i] >max){
				max = fabs(A->entry[j][i]);
			}
		}
	}

	return max;
}

double OrthonormalityCheck(DenseMatrix *Q)
/* Function that calculates the maximum absolute value of entries in
I-Q^T Q. If it returns ~0, the matrix is probably orthonormal. */
{
  //To save memory, I_less_QQt was repeatedly overwritten, rather
  //than creating other matrices, Q^TQ, Q^T or I,
	DenseMatrix I_less_QtQ = CreateNullMatrix(Q->nColumns,Q->nColumns);
  int i;

	QtQ(&I_less_QtQ, Q);

  // We now subtract to get QtQ-I which should be the same (we are finding
  // absolute values anyway).
	for(i = 0; i < I_less_QtQ.nColumns;i++){
		I_less_QtQ.entry[i][i] --;
	}

  //Error is calculated, memory freed and program ends.
  double err = MaxMatValue(&I_less_QtQ);
	FreeMatrix(&I_less_QtQ);
  return err;
}

void QtQ(DenseMatrix *destin, DenseMatrix *Q)
/* Function to find Q^tQ for a matrix Q, to be stored in a DenseMatrix destin*/
{
	int i, j, k;

  // To improve cache performance we calculate the transpose of QtQ, but since
  // QtQ is symmetric we do not need to tranpose it back.
	for ( i = 0; i < destin->nRows ; i++) {
    for ( j = 0; j < destin->nColumns; j++) {
      for ( k = 0; k < Q->nRows; k++) {
          destin->entry[i][j] += Q->entry[i][k]*Q->entry[j][k];
      }
    }
  }
}

double Max_Error(DenseMatrix *A, DenseMatrix *QR)
/* Function to find the maximum absolute value of a A-QR */
{
  double max_error = 0, err;
  int i,j;

  // We iterate through columns first
  for ( j = 0; j < QR->nRows; j++) {
    for ( i = 0; i < QR->nColumns; i++) {
      err = fabs(A->entry[i][j] - QR->entry[i][j]);
      if (err > max_error) {
        max_error = err;
      }
    }
  }
  return max_error;
}
