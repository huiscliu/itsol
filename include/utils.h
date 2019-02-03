#ifndef __ITSOL_UTILS_H__
#define __ITSOL_UTILS_H__

#include <stdarg.h>
#include "type-defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* sets.c */
int nnz_arms(arms PreSt,  FILE *ft);
void itsol_errexit(char *f_str, ...);
void * itsol_malloc(int nbytes, char *msg); 
int setupCS(csptr amat, int len, int job); 
int cleanCS(csptr amat);
int nnz_cs (csptr A) ;
int cscpy(csptr amat, csptr bmat);
int setupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E);
int setupVBMat(vbsptr vbmat, int n, int *nB);
int setupILUT(ilutptr amat, int len);
int cleanVBMat(vbsptr vbmat); 
int nnzVBMat(vbsptr vbmat) ;
int memVBMat(vbsptr vbmat); 
int setupVBILU(vbiluptr lu, int n, int *bsz);
int cleanVBILU(vbiluptr lu); 
int cleanILU(iluptr lu);
int cleanILUT(ilutptr amat, int indic);
int cleanP4(p4ptr amat);
int mallocVBRow(vbiluptr lu, int nrow); 
int mallocRow(iluptr lu, int nrow);
void itsol_zrmC(int m, int n, BData data); 
void itsol_copyBData(int m, int n, BData dst, BData src, int isig); 
int CSRcs(int n, double *a, int *ja, int *ia, csptr mat, int rsa); 
int csrvbsrC(int job, int nBlk, int *nB, csptr csmat, vbsptr vbmat);  
int col2vbcol(int col, vbsptr vbmat);
int nnz_vbilu(vbiluptr lu); 
int nnz_lev4(p4ptr levmat, int *lev, FILE *ft);
int setupILU(iluptr lu, int n);
int CS2lum(int n, csptr Amat, iluptr mat, int typ);
int COOcs(int n, int nnz,  double *a, int *ja, int *ia, csptr bmat);
void coocsr_(int*, int*, double*, int*, int*, double*, int*, int*);

int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F, csptr E, csptr C);
void itsol_setup_arms (arms Levmat);
int cleanARMS(arms ArmsPre);
void itsol_coocsc(int n, int nnz, double *val, int *col, int *row, double **a, int **ja, int **ia, int job);
int nnz_ilu(iluptr lu);
int CSClum(int n, double *a, int *ja, int *ia, iluptr mat, int rsa);
int CSClumC(csptr amat, iluptr mat, int rsa);

/* misc.c */
int SparTran(csptr amat, csptr bmat, int job, int flag); 
int coscalC(csptr mata, double *diag, int nrm);
void itsol_dscale(int n, double *dd, double *x, double * y);
void itsol_hilosort(csptr mat, int abval, int hilo);
void itsol_printmat(FILE *ft, csptr A, int i0, int i1);
void itsol_qqsort(int *ja, double *ma, int left, int right);
void itsol_qsort2C(int *ja, double *ma, int left, int right, int abval); 
void itsol_qsort3i(int *wa, int *cor1, int *cor2, int left, int right); 
void itsol_qsortC(int *ja, double *ma, int left, int right, int abval); 
void itsol_qsortR2I(double *wa, int *cor1, int *cor2, int left, int right); 
int qsplitC(double *a, int *ind, int n, int ncut);
int roscalC(csptr mata, double *diag, int nrm);
void itsol_swapj(int v[], int i, int j);
void itsol_swapm(double v[], int i, int j);
double sys_timer(void);
int dumpArmsMat(arms PreSt, FILE *ft);
int outputLU(iluptr lu, char *filename);
int checkperm(int *p, int n);
void itsol_qsortR1I(double *wa, int *cor1, int left, int right);

#ifdef __cplusplus
}
#endif
#endif 
