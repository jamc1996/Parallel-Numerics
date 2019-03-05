#ifndef QRFACTORING_H
#define QRFACTORING_H

/* Header file for qr factorization methods of qrfactoring.c */

void gs(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R);


void hhorth(DenseMatrix *A, DenseMatrix *Q, DenseMatrix *R);

void Pk_timesx(Vector rk, Vector muk, int n);
void Calculate_Omega(Vector mu_k, Vector A_k, int k, int n);
void Pk_TimesX(Vector rk, Vector muk, int n, int k);


#endif
