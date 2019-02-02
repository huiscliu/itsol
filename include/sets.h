#ifndef __ITSOL_SETS_H__
#define __ITSOL_SETS_H__

#include <stdio.h>
#include "type-defs.h"
#include "ext-protos.h"

/* sets.c */
int nnz_arms(arms PreSt,  FILE *ft);
void errexit(char *f_str, ...);
void *Malloc(int nbytes, char *msg); 
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
int cleanILU( iluptr lu );
int cleanILUT(ilutptr amat, int indic);
int cleanP4(p4ptr amat);
int mallocVBRow(vbiluptr lu, int nrow); 
int mallocRow( iluptr lu, int nrow );
void zrmC(int m, int n, BData data); 
void copyBData(int m, int n, BData dst, BData src, int isig); 
int CSRcs(int n, double *a, int *ja, int *ia, csptr mat, int rsa); 
int csrvbsrC(int job, int nBlk, int *nB, csptr csmat, vbsptr vbmat);  
int col2vbcol( int col, vbsptr vbmat );
int nnz_vbilu(vbiluptr lu); 
int nnz_lev4(p4ptr levmat, int *lev, FILE *ft);
int setupILU( iluptr lu, int n );
int CS2lum( int n, csptr Amat, iluptr mat, int typ);
int COOcs(int n, int nnz,  double *a, int *ja, int *ia, csptr bmat);
void coocsr_(int*, int*, double*, int*, int*, double*, int*, int*);

#endif 
