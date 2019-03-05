#ifndef DENSE_H
#define DENSE_H

/*	dense.h -- header file with column major dense matrix structure
 *
 * 	Author: John Cormican
 */

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
