#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "hhorth.h"




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
 //   if (k==4){break;
//	for (i=4;i<10;i++){printf("%lf\n",R->entry[k][i]);}}
    // Omega calulated to zero subdiagonal entries of k-th column of R.
    Calculate_Omega(Omega.entry[k],R->entry[k],k,R->nRows);

    // R_k updated according to P_k R_k
    Pk_TimesX(R->entry[k],Omega.entry[k],R->nRows,k);

    //Q columns calculated
    
  }

}

DenseMatrix CreateLT(int n, int m)
/* Creates enough room to efficiently store non zero values of m vectors with
zeros above the diagonal. */
{
  // Storage space calculated

  DenseMatrix A;
	int sum = m*(n-m);
  int i;
  for ( i = 0; i < m; i++) {
    sum += i;
  }
	A.nRows = n;
	A.nColumns = m;
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

void mod_hhorth(DenseMatrix *A, DenseMatrix *W, DenseMatrix *R)
/* Implements QR factorization of matrix A via Householder orthogonalisation.
Updates the values of dense matrices Q and R. */
{
  //Iterators declared and place to store omega
  int i, k;

  //DenseMatrix Omega = CreateNullMatrix(Q->nRows, Q->nColumns);


  for ( k = 0; k < R->nColumns; k++)
  {
    // Values of A_k are copied into R_k
    memcpy(R->entry[k],A->entry[k],(R->nRows)*sizeof(double));
    // Loop takes into account the effects of previous P_k transformations.
    if (k>0)
    {
        for (i = 0; i < k; i++)
        {
          Pk_TimesX(R->entry[k],W->entry[i],R->nRows,i);
        }
    }

    // Omega calulated to zero k-th column of R.
    Calculate_Omega(W->entry[k],R->entry[k],k,R->nRows);


    Pk_TimesX(R->entry[k],W->entry[k],R->nRows,k);
  }

}

void mod_mod_hhorth(DenseMatrix *A, DenseMatrix *W, DenseMatrix *R, int col)
/* Implements QR factorization of matrix A via Householder orthogonalisation.
Updates the values of dense matrices Q and R. */
{
  //Iterators declared and place to store omega
  int i, k;

  //DenseMatrix Omega = CreateNullMatrix(Q->nRows, Q->nColumns);



  for ( k = 0; k < R->nColumns; k++)
  {
    // Values of A_k are copied into R_k
    memcpy(R->entry[k],A->entry[k],(R->nRows)*sizeof(double));
    // Loop takes into account the effects of previous P_k transformations.
    if (k>0)
    {
        for (i = 0; i < k; i++)
        {
          Pk_TimesX(R->entry[k],W->entry[i],R->nRows,i);
        }
    }

    // Omega calulated to zero k-th column of R.
    Calculate_Omega(W->entry[k],R->entry[k],k,R->nRows);


    Pk_TimesX(R->entry[k],W->entry[k],R->nRows,k);
  }

}

void ApplyQT(DenseMatrix* B,DenseMatrix* Q,DenseMatrix* W)
{
	int i,j;
	memcpy(B->data_,W->data_,sizeof(double)*B->nRows*B->nColumns);
	for (i=0; i<B->nColumns;i++)
	{
		for (j=0; j<Q->nColumns; j++)
		{
			Pk_TimesX(B->entry[i],Q->entry[j],B->nRows,j);
		}
	}
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

double two_norm(Vector Aj,int m){
  int i;
  double sum = 0;
  for ( i = 0; i < m; i++) {
    sum+= Aj[i]*Aj[i];
  }
  return sqrt(sum);
}

DenseMatrix CreateNullMatrix(int n, int m)
{
  DenseMatrix A;
  int j;
  A.nRows = n;
  A.nColumns = m;
  A.data_ = malloc(sizeof(double)*n*m);//all entries stored in t.data_.
  memset(A.data_,0.0,sizeof(double)*n*m);
  A.entry = malloc(sizeof(double*)*m);//space allocated for m columns.
  for (j=0;j<m;j++)
  {
    A.entry[j] = &A.data_[j*n];
  }
  return A;
}

void SetToRandom(DenseMatrix *A)
{
  int i;
  for (i=0;i<A->nRows*A->nColumns ;i++)
  {
    A->data_[i]= (drand48()*4) - 2;
  }
}

void FreeMatrix(DenseMatrix *A)
{
  free(A->data_);
  free(A->entry);
}


void PrintMatrix(DenseMatrix *A){
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
}
