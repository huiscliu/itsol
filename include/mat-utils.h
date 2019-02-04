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
int itsol_invGauss(int nn, double *A);

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
int itsol_invSVD(int nn, double *A);

void itsol_matvecC(ITS_CsPtr mat, double *x, double *y);
void itsol_matvecCSC(SMatptr mat, double *x, double *y);
int itsol_CondestC(ITS_IluPtr lu, FILE * fp);
int itsol_diag_scal(ITS_VbsPtr vbmat);
int itsol_diagvec(ITS_VbsPtr vbmat, BData x, BData y);
void itsol_matvec(ITS_CsPtr mata, double *x, double *y); 
void itsol_matvecCSR(SMatptr mat, double *x, double *y);
void itsol_matvecz(ITS_CsPtr mata, double *x, double *y, double *z);
void itsol_vbmatvec(ITS_VbsPtr vbmat, double *x, double *y);
void itsol_luinv(int n, double *a, double *x, double *y); 
int itsol_vblusolC(double *y, double *x, ITS_VbiluPtr lu); 
int itsol_lusolC(double *y, double *x, ITS_IluPtr lu); 
int itsol_rpermC(ITS_CsPtr mat, int *perm); 
int itsol_cpermC(ITS_CsPtr mat, int *perm) ; 
int itsol_dpermC(ITS_CsPtr mat, int *perm) ; 
int itsol_CSparTran(ITS_CsPtr amat, ITS_CsPtr bmat, CompressType *compress);
double itsol_vbnorm2(int sz, double *a);
void itsol_Lsol(ITS_CsPtr mata, double *b, double *x);
void itsol_Usol(ITS_CsPtr mata, double *b, double *x);
int itsol_ascend (ITS_P4Ptr levmat, double *x, double *wk);
int itsol_descend(ITS_P4Ptr levmat, double *x, double *wk);
int itsol_armsol2(double *x, arms Prec);
int itsol_condestArms(arms armspre, double *y, FILE *fp);
int itsol_VBcondestC(ITS_VbiluPtr, FILE *fp); 
int itsol_CondestLUM(ITS_IluPtr lu, double *y, double *x, FILE *fp);
void itsol_matvecVBR(SMatptr mat, double *x, double *y);
void itsol_matvecLDU(SMatptr mat, double *x, double *y);
int itsol_preconILU(double *x, double *y, SPreptr mat);
int itsol_preconVBR(double *x, double *y, SPreptr mat);
int itsol_preconLDU(double *x, double *y, SPreptr mat);
int itsol_preconARMS(double *x, double *y, SPreptr mat);
ITS_P4Ptr itsol_Lvsol2(double *x, int nlev, ITS_P4Ptr levmat, ilutptr ilusch) ;
int itsol_Uvsol2(double *x, int nlev, int n, ITS_P4Ptr levmat, ilutptr ilusch); 
void itsol_SchLsol(ilutptr ilusch, double *y) ;
void itsol_SchUsol(ilutptr ilusch, double *y) ;
int itsol_lumsolC(double *y, double *x, ITS_IluPtr lu);
int itsol_condestLU(ITS_IluPtr lu, FILE *fp);
int itsol_invGauss(int nn, double *A); 

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
int itsol_init_blocks(ITS_CsPtr csmat, int *pnBlock, int **pnB, int **pperm,
        double eps, double *t_hash, double *t_angle);

#ifdef __cplusplus
}
#endif

#endif 
