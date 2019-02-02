#include "type-defs.h"

/* *-------------------- inversion by svd
   This calls lapack routines for inverting a dense matrix.
   dgetrf and dgetri

   ON ENTRY
 ** A = square matrix of size n x n -- dimensioned with
 ** leading dimension =  n

 ON RETURN A contains the inverse  of the input matrix.
 */
int invGauss(int nn, double *A);

/* *-------------------- inversion by svd
   This calls lapack routine dgesvd --
   ON ENTRY
 ** A = square matrix of size n x n -- dimensioned with
 ** leading dimension =  n
 ON RETURN A contains the truncated SVD inverse of input matrix.
 ** tolerance set for truncation is TOL and can be changed in
 ** above define statement
 **--------------------
 */
int invSVD(int nn, double *A);
