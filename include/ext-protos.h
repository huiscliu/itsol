#ifndef __ITSOL_EXT_PROTOS_H__
#define __ITSOL_EXT_PROTOS_H__

#include <stdio.h>
#include "type-defs.h"

#if defined(_SGI) || defined(_LINUX)
#define dnrm2   dnrm2_
#define ddot    ddot_
#define daxpy   daxpy_
#define qsplit  qsplit_
#define dscal   dscal_
#define dgemv   dgemv_
#define dgemm   dgemm_
#define dgetrf  dgetrf_
#define dgetri  dgetri_
#define dgesvd  dgesvd_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#define qsplit  qsplit_
#elif defined(_IBM)
#include <essl.h>
#define dnrm2   dnrm2
#define ddot    ddot
#define daxpy   daxpy
#define qsplit  qsplit
#define dscal   dscal
#define dgemv   dgemv
#define dgemm   dgemm
#define dgetrf  dgetrf
#define dgetri  dgetri
#define dgesvd  dgesvd
#define readmtc readmtc
#define csrcsc  csrcsc 
#define roscal  roscal
#define coscal  coscal
#define qsplit  qsplit
#else
#define dnrm2   dnrm2_
#define ddot    ddot_
#define daxpy   daxpy_
#define qsplit  qsplit_
#define dscal   dscal_
#define dgemv   dgemv_
#define dgemm   dgemm_
#define dgetrf  dgetrf_
#define dgetri  dgetri_
#define dgesvd  dgesvd_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#define qsplit  qsplit_
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

/* FORTRAN routines */
void readmtc(int*,  int*,  int*,  char*,  double*,  int*,
	     int*,  double*, int*,  char*,  int*,  int*,  int*,
	     char*,  char*, char*,  int*) ;
void csrcsc(int*, int*, int*, double*, int*, int*, double*,
		    int*, int*) ; 
void qsplit(double *a, int *ind, int *n, int *ncut);	
void dgesvd(char*, char*, int*, int*, double*, int*, double*,
		   double *, int*, double*, int*, double*, int*,
		   int*); 
void csrcoo( int *, int *, int *, double *, int *, int *, int *,
		double *, int *, int *, int *);    

#if defined(_IBM)
#define DDOT(n,x,incx,y,incy)        ddot((n), (x), (incx), (y), (incy)) 
#define DCOPY(n,x,incx,y,incy)       dcopy((n), (x), (incx), (y), \
					   (incy)) 
#define DSCAL(n,alpha,x,incx)        dscal((n), (alpha), (x), (incx)) 
#define DAXPY(n,alpha,x,incx,y,incy) daxpy((n), (alpha), (x), (incx), (y), (incy)) 
#define DNRM2(n,x,incx)              dnrm2((n), (x), (incx))

#define IDMIN(n,sx,incx)             idmin((n), (sx), (incx))
#define DGEMV(transa,m,n,alpha,a,lda,x,incx,beta,y,incy)		\
  dgemv((transa), (m), (n),						\
	(alpha), (a), (lda), (x), (incx),				\
	(beta), (y), (incy))

#define DGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)		\
  dgemm((transa),(transb),						\
	(l),(n),(m),(alpha),(a),					\
	(lda),(b),(ldb),(beta),(c),(ldc))
#define DGETRF(m, n, a, lda, ipvt, info)  \
  dgetrf((m), (n), (a), (lda), (ipvt), (info))
#define DGETRI(n, a, lda, ipvt, work, lwork, info)		\
  dgetri((n), (a), (lda), (ipvt), (work), (lwork), (info))
#else
#define DDOT(n,x,incx,y,incy)        ddot(&(n),(x),&(incx),(y),&(incy))
#define DCOPY(n,x,incx,y,incy)       dcopy(&(n),(x),&(incx),(y),&(incy))
#define DSCAL(n,alpha,x,incx)        dscal(&(n),&(alpha),(x), &(incx))
#define DAXPY(n,alpha,x,incx,y,incy) daxpy(&(n), &(alpha), (x), &(incx), y, &(incy))
#define DNRM2(n, x, incx)            dnrm2(&(n), (x), &(incx))
#define IDMIN(n, sx, incx)           idmin((&(n), (sx), &(incx))
#define DGEMV(transa, m, n, alpha, a, lda, x, incx, beta, y, incy)  \
  dgemv((transa), &(m), &(n), &(alpha), (a), &(lda), (x), &(incx), \
	 &(beta), (y), &(incy))
#define DGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)	\
  dgemm((transa), (transb), &(l), &(n), &(m), &(alpha), (a),	\
	 &(lda), b, &(ldb), &(beta), (c), &(ldc)) 
#define DGETRF(m, n, a, lda, ipvt, info)		\
  dgetrf(&(m), &(n), (a), &(lda), (ipvt), (info))
#define DGETRI(n, a, lda, ipvt, work, lwork, info)			\
  dgetri(&(n), (a), &(lda), (ipvt), (work), &(lwork), (info))

double ddot(int *n, double *x, int *incx, double *y, int *incy);  
void   dcopy(int *n, double *x, int *incx, double *y, int *incy); 
void   dscal(int *n, double *alpha, double *x, int *incx);
void   daxpy(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
double dnrm2(int *n, double *x, int *incx);
void   idmin(int *n, double *sx, int *incx);
void   dgemv(char *transa, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx,
    double *beta, double *y, int *incy);

void   dgemm(char *transa, char *transb, int *l, int *m, int *n, double *alpha, double *a, int *lda,
        double *b, int *ldb, double *beta, double *c, int *ldc);       

void   dgetrf(int *m, int *n, double *a, int *lda, int *ipvt, int *info); 
void   dgetri(int *n, double *a, int *lda, int *ipvt, double *work, int *lwork, int *info);
#endif 

#endif 
