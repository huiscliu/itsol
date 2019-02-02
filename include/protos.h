#ifndef __ITSOL_INCLUDED_PROTOS_H__
#define __ITSOL_INCLUDED_PROTOS_H__

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
int invSVD(int nn, double *A) ;

/* setblks.c */ 
int KeyComp(const void *vfst, const void *vsnd);
int init_blocks(csptr, int *, int **, int **, double, double *, double *);  

/* upper directory */
int vbilukC( int lofM, vbsptr vbmat, vbiluptr lu, FILE *fp ); 
int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol, BData *w, FILE *fp ); 
int ilutc(iluptr mt, iluptr lu, int lfil, double tol, int drop, FILE *fp ); 
int ilukC( int lofM, csptr csmat, iluptr lu, FILE *fp );
int ilut( csptr csmat, iluptr lu, int lfil, double tol, FILE *fp );
int fgmr(SMatptr Amat, SPreptr lu, double *rhs, double *sol, double tol,
        int im, int *itmax, FILE *fits ); 
int arms2(csptr Amat, int *ipar, double *droptol, int *lfil, 
        double tolind, arms PreMat, FILE *ft) ;
int condestLU( iluptr, FILE *);
int nnz_ilu( iluptr lu ); 
void roscal (int* nrow, int* job, int* nrm, double *a, int *ja,
        int *ia, double *diag, double *b, int *jb, int
        *ib, int *ierr) ;  
void coscal (int* nrow, int* job, int* nrm, double *a, int *ja,
        int *ia, double *diag, double *b, int *jb, int *ib, int *ierr) ;  
int outputLU( iluptr lu, char *filename );
int lumsolC(double *y, double *x, iluptr lu );
void lumatvec(iluptr mat, double *x, double *y);
int CSClum( int n, double *a, int *ja, int *ia, iluptr mat, int rsa ); 
int CSClumC(csptr amat, iluptr mat, int rsa ); 
void setup_arms (arms Levmat);
int cleanARMS(arms ArmsPre);
int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F, csptr E, csptr C); 

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

/* piluNEW.c */
int pilu(p4ptr amat, csptr B, csptr C, double *droptol, int *lfil, csptr schur);

/* ilutpC.c */
int ilutD(csptr amat, double *droptol, int *lfil, ilutptr ilusch);
int ilutpC(csptr amat, double *droptol, int *lfil, double
        permtol, int mband, ilutptr ilusch);

/* PQ.c */
int PQperm(csptr mat, int bsize, int *Pord, int *Qord, int *nnod, double tol);
int add2com(int *nback, int nod, int *iord, int *riord);
int add2is(int *last, int nod, int *iord, int *riord);
int indsetC(csptr mat, int bsize, int *iord, int *nnod, double tol); 
int preSel(csptr mat, int *icor, int *jcor, int job, double tol, int *count);

/* indsetC.c */
int weightsC(csptr mat, double *w);

/* setblks.c */
int KeyComp( const void *vfst, const void *vsnd );
int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm, double eps, double *t_hash,
        double *t_angle );

/* systimer.c */
double sys_timer(void);

/*auxill.c */
int read_inputs( char *in_file, io_t *pio );
int get_matrix_info( FILE *fmat, io_t *pio );
void output_blocks( int nBlock, int *nB, FILE *f );
void output_perm( int n, int *perm, FILE *f );
int read_coo(double **VAL, int **COL, int **ROW, io_t *pio, double **rhs, double **sol, int job);
int readhb_c(int *NN, double **AA, int **JA, int **IA, io_t *pio, double **rhs, double **sol, int *rsa);
int readhb_2(int *NN, double **AA, int **JA, int **IA, io_t *pio, double **rhs, double **sol, int *rsa, int fmt);
void output_header( io_t *pio );
void output_header_vb( io_t *pio );
void output_result( int lfil, io_t *pio, int iparam );
void set_arms_pars(io_t* io, int Dscale, int *ipar, double *dropcoef, int *lfil);
void randvec (double *v, int n);

int lutsolC( double *y, double *x, iluptr lu );

int dumpArmsMat(arms PreSt, FILE *ft);
int checkperm(int *p, int n) ;
void qsortR1I(double *wa, int *cor1, int left, int right);
void swapj(int v[], int i, int j);
void swapm(double v[], int i, int j) ;
void coocsc(int n, int nnz, double *val, int *col, int *row, double **a, int **ja, int **ia, int job);

#endif 
