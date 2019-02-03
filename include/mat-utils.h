#ifndef __ITSOL_MatOps_H__
#define __ITSOL_MatOps_H__

#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

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

void matvecC(csptr mat, double *x, double *y);
void matvecCSC(SMatptr mat, double *x, double *y);
int CondestC(iluptr lu, FILE * fp);
int diag_scal(vbsptr vbmat);
int diagvec(vbsptr vbmat, BData x, BData y);
void matvec(csptr mata, double *x, double *y); 
void matvecCSR(SMatptr mat, double *x, double *y);
void matvecz(csptr mata, double *x, double *y, double *z);
void vbmatvec(vbsptr vbmat, double *x, double *y);
void luinv(int n, double *a, double *x, double *y); 
int vblusolC(double *y, double *x, vbiluptr lu); 
int lusolC(double *y, double *x, iluptr lu); 
int rpermC(csptr mat, int *perm); 
int cpermC(csptr mat, int *perm) ; 
int dpermC(csptr mat, int *perm) ; 
int CSparTran(csptr amat, csptr bmat, CompressType *compress);
double vbnorm2(int sz, double *a);
void Lsol(csptr mata, double *b, double *x);
void Usol(csptr mata, double *b, double *x);
int ascend (p4ptr levmat, double *x, double *wk);
int descend(p4ptr levmat, double *x, double *wk);
int armsol2(double *x, arms Prec);
int condestArms(arms armspre, double *y, FILE *fp);
int VBcondestC(vbiluptr, FILE *fp); 
int CondestLUM(iluptr lu, double *y, double *x, FILE *fp);
void matvecVBR(SMatptr mat, double *x, double *y);
void matvecLDU(SMatptr mat, double *x, double *y);
int preconILU(double *x, double *y, SPreptr mat);
int preconVBR(double *x, double *y, SPreptr mat);
int preconLDU(double *x, double *y, SPreptr mat);
int preconARMS(double *x, double *y, SPreptr mat);
p4ptr itsol_Lvsol2(double *x, int nlev, p4ptr levmat, ilutptr ilusch) ;
int Uvsol2(double *x, int nlev, int n, p4ptr levmat, ilutptr ilusch); 
void SchLsol(ilutptr ilusch, double *y) ;
void SchUsol(ilutptr ilusch, double *y) ;
int lumsolC(double *y, double *x, iluptr lu);
int condestLU(iluptr lu, FILE *fp);
int invGauss(int nn, double *A); 

/*----------------------------------------------------------------------------
 * Setup Blocks (rows and columns might be permuted to get better results)
 *----------------------------------------------------------------------------
 * Na Li, Aug 2001
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * csmat   = a matrix stored in SpaFmt format
 * eps     = parameter for deciding when to do a union of two rows
 *           into the same group.  Two rows u and v are merged into a 
 *           block  when cos(<u,v>) == (u,v)/(|u|*|v|), is > eps. 
 *           eps should be <= 1. 
 *----------------------------------------------------------------------------
 * on return:
 * ==========
 * csmat   = matrix stored in SpaFmt format after permutation
 * pnBlock = dimension of the block matrix
 * pnB     = dimension of each block
 *
 *----------------------------------------------------------------------------
 * Combination of hash method and angle method:
 *----------------------------------------------------------------------------
 * Designed for the matrices with symmetric patterns
 * (1) Hash method
 *     a. Calculate hash values
 *     b. qsort rows according to their hash values
 *     c. Get compressed graph as the following format:
 * (2) Angle method
 *     a. Calculate A^T
 *     b. for i-th row, calculate dot product (row_i, row_j) using A*A^T
 *        algorithm where j = i+1, ..., n-1 and group[j] == -1
 *        if cos(<row_i, row_j>) = (row_i,row_j)/|row_i||row_j| is > eps,
 *        we merge row_i and row_j by resetting
 *        group[j] = i and size[i] = size[i]+size[j]
 *--------------------------------------------------------------------------*/
int itsol_init_blocks(csptr csmat, int *pnBlock, int **pnB, int **pperm,
        double eps, double *t_hash, double *t_angle);

#ifdef __cplusplus
}
#endif

#endif 
