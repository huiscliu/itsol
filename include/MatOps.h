#ifndef __ITSOL_MatOps_H__
#define __ITSOL_MatOps_H__

#include "sets.h"
#include "svdInvC.h"

/* MatOps.c */
void matvecC(csptr mat, double *x, double *y );
void matvecCSC(SMatptr mat, double *x, double *y);
int CondestC(iluptr lu, FILE *fp );
int diag_scal(vbsptr vbmat);
int diagvec(vbsptr vbmat, BData x, BData y);
void matvec(csptr mata, double *x, double *y); 
void matvecCSR(SMatptr mat, double *x, double *y);
void matvecz(csptr mata, double *x, double *y, double *z);
void vbmatvec(vbsptr vbmat, double *x, double *y);
void luinv(int n, double *a, double *x, double *y); 
int vblusolC(double *y, double *x, vbiluptr lu); 
int lusolC( double *y, double *x, iluptr lu ); 
int rpermC(csptr mat, int *perm); 
int cpermC(csptr mat, int *perm) ; 
int dpermC(csptr mat, int *perm) ; 
int CSparTran(csptr amat, csptr bmat, CompressType *compress);
double vbnorm2(int sz, double *a);
void Lsol(csptr mata, double *b, double *x);
void Usol(csptr mata, double *b, double *x);
int ascend (p4ptr levmat, double *x, double *wk);
int descend(p4ptr levmat, double *x, double *wk);
int armsol2(double *x,  arms Prec);
int condestArms(arms armspre, double *y, FILE *fp );
int VBcondestC( vbiluptr, FILE *fp ); 
int CondestLUM(iluptr lu, double *y, double *x, FILE *fp );
void matvecVBR(SMatptr mat, double *x, double *y);
void matvecLDU(SMatptr mat, double *x, double *y);
int  preconILU(double *x, double *y, SPreptr mat);
int  preconVBR(double *x, double *y, SPreptr mat);
int  preconLDU(double *x, double *y, SPreptr mat);
int  preconARMS(double *x, double *y, SPreptr mat);
p4ptr Lvsol2(double *x, int nlev, p4ptr levmat, ilutptr ilusch) ;
int   Uvsol2(double *x, int nlev, int n, p4ptr levmat, ilutptr
		    ilusch); 
void  SchLsol(ilutptr ilusch, double *y) ;
void  SchUsol(ilutptr ilusch, double *y) ;
int lumsolC(double *y, double *x, iluptr lu );
int condestLU( iluptr lu, FILE *fp );

#ifndef _IBM
void dgemv(char*, int *, int*, double*, double *, int*,
		  double*, int*, double*, double*, int*);  
void dgemm(char*, char*, int*, int*, int*, double*, double*,
		  int*, double*, int*, double*, double*, int*) ;  
void dgetrf(int*, int*, double*, int*, int*, int*); 
void dgetri(int*, double*, int*, int*, double*,  int*, int* );
double dnrm2( int *, double *, int * );
void dscal(int*, double*, double*, int*); 
#endif 

int invGauss(int nn, double *A); 

#endif 
