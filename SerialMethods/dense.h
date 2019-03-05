#ifndef DENSE_H
#define DENSE_H

// Typecast used increase readability
typedef double* Vector;

// DenseMatrix structure organised to ensure contiguous memory.
typedef struct
{
  double*   data_;
  int       nRows;
  int       nColumns;
  double**  entry;
}
DenseMatrix;

#endif
