#ifndef MATOPS_H
#define MATOPS_H

#include "dense.h"

// Operations used in orthogonalisation process
double two_norm(Vector Aj, int m);
void scalar_div(Vector Ai,double Rii, Vector Qj,int j);
double inner_prod(Vector Aj, Vector Qi, int n);
void scalar_mult(Vector Ai,double Rii, Vector Qj,int j);

// Operations used for error checking
void MatrixMultiply(DenseMatrix *Q, DenseMatrix *R, DenseMatrix *QR);
double MaxMatValue(DenseMatrix *A);
double Max_Error(DenseMatrix *A, DenseMatrix *QR);
double OrthonormalityCheck(DenseMatrix *Q);
void QtQ(DenseMatrix *destin, DenseMatrix *Q);

#endif
