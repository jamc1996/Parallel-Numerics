#ifndef HHORTH_H
#define HHORTH_H

/*	hhorth.h -- header file for hhorth.c 
 *
 *	Author: John Cormican
 *
 * */


#include "dense.h"

#include <stdio.h>

void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R);
void mod_hhorth(DenseMatrix *A, DenseMatrix *W, DenseMatrix *R);
void ApplyQT(DenseMatrix* B,DenseMatrix* Q,DenseMatrix* W);
void Pk_timesx(Vector rk, Vector muk, int n);
void Calculate_Omega(Vector mu_k, Vector A_k, int k, int n);
void Pk_TimesX(Vector rk, Vector muk, int n, int k);
DenseMatrix CreateLT(int n, int m);
double two_norm(Vector Aj, int m);

DenseMatrix CreateNullMatrix(int n, int m);
void SetToRandom(DenseMatrix *A);
void FreeMatrix(DenseMatrix *A);
void PrintMatrix(DenseMatrix *A);


#endif
