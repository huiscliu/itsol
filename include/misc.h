#ifndef __ITSOL_MISC_H__
#define __ITSOL_MISC_H__

#include "type-defs.h"

/* misc.c */
int SparTran(csptr amat, csptr bmat, int job, int flag); 
int coscalC(csptr mata, double *diag, int nrm);
void dscale(int n, double *dd, double *x, double * y);
void hilosort(csptr mat, int abval, int hilo);
void printmat(FILE *ft, csptr A, int i0, int i1);
void qqsort(int *ja, double *ma, int left, int right);
void qsort2C(int *ja, double *ma, int left, int right, int abval); 
void qsort3i(int *wa, int *cor1, int *cor2, int left, int right); 
void qsortC(int *ja, double *ma, int left, int right, int abval); 
void qsortR2I(double *wa, int *cor1, int *cor2, int left, int right); 
int qsplitC(double *a, int *ind, int n, int ncut);
int roscalC(csptr mata, double *diag, int nrm);
void swapj(int v[], int i, int j);
void swapm(double v[], int i, int j);
double sys_timer(void);
int dumpArmsMat(arms PreSt, FILE *ft);
int outputLU( iluptr lu, char *filename);
int checkperm(int *p, int n);
void qsortR1I(double *wa, int *cor1, int left, int right);

#endif 
