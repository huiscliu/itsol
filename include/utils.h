#ifndef __ITSOL_UTILS_H__
#define __ITSOL_UTILS_H__

#include <stdarg.h>
#include "type-defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* sets.c */
int itsol_nnz_arms(arms PreSt,  FILE *ft);
void itsol_errexit(char *f_str, ...);
void * itsol_malloc(int nbytes, char *msg); 
int itsol_setupCS(ITS_CsPtr amat, int len, int job); 
int itsol_cleanCS(ITS_CsPtr amat);
int itsol_nnz_cs (ITS_CsPtr A) ;
int itsol_cscpy(ITS_CsPtr amat, ITS_CsPtr bmat);
int itsol_setupP4 (p4ptr amat, int Bn, int Cn,  ITS_CsPtr F,  ITS_CsPtr E);
int itsol_setupVBMat(vbsptr vbmat, int n, int *nB);
int itsol_setupILUT(ilutptr amat, int len);
int itsol_cleanVBMat(vbsptr vbmat); 
int itsol_nnzVBMat(vbsptr vbmat) ;
int itsol_memVBMat(vbsptr vbmat); 
int itsol_setupVBILU(vbiluptr lu, int n, int *bsz);
int itsol_cleanVBILU(vbiluptr lu); 
int itsol_cleanILU(iluptr lu);
int itsol_cleanILUT(ilutptr amat, int indic);
int itsol_cleanP4(p4ptr amat);
int itsol_mallocVBRow(vbiluptr lu, int nrow); 
int itsol_mallocRow(iluptr lu, int nrow);
void itsol_zrmC(int m, int n, BData data); 
void itsol_copyBData(int m, int n, BData dst, BData src, int isig); 
int itsol_CSRcs(int n, double *a, int *ja, int *ia, ITS_CsPtr mat, int rsa); 
int itsol_csrvbsrC(int job, int nBlk, int *nB, ITS_CsPtr csmat, vbsptr vbmat);  
int itsol_col2vbcol(int col, vbsptr vbmat);
int itsol_nnz_vbilu(vbiluptr lu); 
int itsol_nnz_lev4(p4ptr levmat, int *lev, FILE *ft);
int itsol_setupILU(iluptr lu, int n);
int itsol_CS2lum(int n, ITS_CsPtr Amat, iluptr mat, int typ);
int itsol_COOcs(int n, int nnz,  double *a, int *ja, int *ia, ITS_CsPtr bmat);
void itsol_coocsr_(int*, int*, double*, int*, int*, double*, int*, int*);

int itsol_csSplit4(ITS_CsPtr amat, int bsize, int csize, ITS_CsPtr B, ITS_CsPtr F, ITS_CsPtr E, ITS_CsPtr C);
void itsol_setup_arms (arms Levmat);
int itsol_cleanARMS(arms ArmsPre);
void itsol_coocsc(int n, int nnz, double *val, int *col, int *row, double **a, int **ja, int **ia, int job);
int itsol_nnz_ilu(iluptr lu);
int itsol_CSClum(int n, double *a, int *ja, int *ia, iluptr mat, int rsa);
int itsol_CSClumC(ITS_CsPtr amat, iluptr mat, int rsa);

/* misc.c */
int itsol_SparTran(ITS_CsPtr amat, ITS_CsPtr bmat, int job, int flag); 
int itsol_coscalC(ITS_CsPtr mata, double *diag, int nrm);
void itsol_dscale(int n, double *dd, double *x, double * y);
void itsol_hilosort(ITS_CsPtr mat, int abval, int hilo);
void itsol_printmat(FILE *ft, ITS_CsPtr A, int i0, int i1);
void itsol_qqsort(int *ja, double *ma, int left, int right);
void itsol_qsort2C(int *ja, double *ma, int left, int right, int abval); 
void itsol_qsort3i(int *wa, int *cor1, int *cor2, int left, int right); 
void itsol_qsortC(int *ja, double *ma, int left, int right, int abval); 
void itsol_qsortR2I(double *wa, int *cor1, int *cor2, int left, int right); 
int itsol_qsplitC(double *a, int *ind, int n, int ncut);
int itsol_roscalC(ITS_CsPtr mata, double *diag, int nrm);
void itsol_swapj(int v[], int i, int j);
void itsol_swapm(double v[], int i, int j);
double itsol_get_time(void);
int itsol_dumpArmsMat(arms PreSt, FILE *ft);
int itsol_outputLU(iluptr lu, char *filename);
int itsol_checkperm(int *p, int n);
void itsol_qsortR1I(double *wa, int *cor1, int left, int right);

#ifdef __cplusplus
}
#endif
#endif 
