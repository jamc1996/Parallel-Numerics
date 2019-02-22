#ifndef DENSE_H
#define DENSE_H

typedef double* Vector;
typedef struct
{
  double*   data_;
  int       nRows;
  int       nColumns;
  double**  entry;
}
DenseMatrix;

#endif
