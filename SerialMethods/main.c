#include <stdio.h>

#include "matops.h"
#include "matSetUp.h"
#include "dense.h"
#include "qrfactoring.h"

/*
		main.c -- 	program to test the functions for qr factorisation found in the
                qrfactoring.c using the functions of matop.c .

		Author:			John Cormican

		Purpose:		To test the functions provided in qrfactoring.c.

		Usage:			Run main. main calls run_test() calls the functions for qr
                factorisation and writes the results to files.
*/


int run_test(int n, int m, int num_tests, char *filename, char* filename2);

int main(int argc, char *argv[])
/* Begins program by testing results of qr factorisation by calling run_test*/
{
  int n = 100;
  int m = 90;
  int num_tests = 10;

  run_test(n,m,num_tests,"MatrixErr.txt","Orthonorm.txt");

  return 0;
}

int run_test(int n, int m, int num_tests, char *filename, char *filename2)
/* Function creates nx, matrices for num_tests number of iterations, testing
the success of qr factorisation and printing the results to files. */
{
  // Variables and Matrices Allocated.
  FILE *fp, *fp2;
  fp = fopen(filename,"w");
  fp2 = fopen(filename2,"w");
  int i;
  double err;
  DenseMatrix A = CreateRandomMatrix(n,m);
  DenseMatrix Q = CreateNullMatrix(n,m);
  DenseMatrix R = CreateNullMatrix(n,m);
  DenseMatrix QR = CreateNullMatrix(n,m);

  fprintf(fp,"|| Test No.  || Error (A-QR) using basic gs || Error (A-QR) using householder ||\n" );
  fprintf(fp2,"|| Test No.  || I - (Q^T)Q using basic gs || I - (Q^T)Q using householder ||\n" );

  // Now qr factorisation performed num_tests times
  for ( i = 0; i<num_tests; i++)
  {
    fprintf(fp, "||    %02d     ||",i+1 );
    fprintf(fp2, "||    %02d     ||",i+1 );


    // A entries set to new random values, Q,R,QR set to zeros.
	  SetToRandom(&A);
    SetToZeros(&QR);
    SetToZeros(&Q);
    SetToZeros(&R);

    gs(&A,&Q,&R);

    // Now test if A = QR and that Q is orthonormal.

    MatrixMultiply(&Q,&R,&QR);
    err = Max_Error(&A,&QR);
    fprintf(fp,"        %E         ||",err );

		err = OrthonormalityCheck(&Q);
    fprintf(fp2,"       %E        ||",err );



    // Same procedure performed on A but using householder orthogonalisation.

    SetToZeros(&QR);
    SetToZeros(&Q);
    SetToZeros(&R);

    hhorth(&A,&Q,&R);

    MatrixMultiply(&Q,&R,&QR);
    err = 0;
    err = Max_Error(&A,&QR);
    fprintf(fp,"          %E          ||\n",err );

    err = OrthonormalityCheck(&Q);
    fprintf(fp2,"         %E         ||\n",err );

  }
  printf("Results of tests printed to files %s and %s.\n",filename,filename2 );

  // Memory freed, files closed and program ends
	FreeMatrix(&A);
	FreeMatrix(&Q);
	FreeMatrix(&R);
  FreeMatrix(&QR);
  fclose(fp);
  fclose(fp2);
  return 0;
}
