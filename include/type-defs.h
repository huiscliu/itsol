#ifndef __VBLOCK_TYPE_DEFS_H__
#define __VBLOCK_TYPE_DEFS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h> 

#include "config.h"
#include "ext-protos.h"

#define ITS_MAX_BLOCK_SIZE   100

/* FORTRAN style vblock format, compatible for many FORTRAN routines */
#define ITS_DATA(a,row,i,j)  (a[(j)*(row)+(i)])

/* the dimension of ith Block */
#define ITS_B_DIM(bs,i)      (bs[i+1]-bs[i])

#define ITS_MAX_LINE        256
#define ITS_MAX_HBNAME      64
#define ITS_MAX_MAT	        100
#define ITS_MaxNamLen       64
#define ITS_HB   1
#define ITS_MM0  2
#define ITS_MM1  3
#define ITS_UNK  4

/*--------------------------------------------- 
  | C-style CSR format - used internally
  | for all matrices in CSR format 
  |---------------------------------------------*/
typedef struct ITS_SparMat_
{
    int n;
    int *nzcount;  /* length of each row */
    int **ja;      /* pointer-to-pointer to store column indices  */
    double **ma;   /* pointer-to-pointer to store nonzero entries */

} ITS_SparMat, *ITS_CsPtr;

typedef double *ITS_BData;

typedef struct ITS_VBSparMat_
{
    int n;        /* the block row dimension of the matrix      */
    int *bsz;     /* the row/col of the first element of each   */

    /* diagonal block                             */
    int *nzcount;     /* length of each row                         */
    int **ja;         /* pointer-to-pointer to store column indices */
    ITS_BData **ba;   /* pointer-to-pointer to store nonzero blocks */
    ITS_BData *D;     /* to store inversion of diagonals            */

} ITS_VBSparMat, *ITS_VBSPtr;

typedef struct ITS_VBILUSpar_
{
    int n;
    int *bsz;     /* the row/col of the first element of each   */

    /* diagonal block                             */
    ITS_BData *D;     /* diagonal blocks                            */
    ITS_VBSPtr L;     /* L part blocks                              */
    ITS_VBSPtr U;     /* U part blocks                              */
    int *work;    /* working buffer                             */
    ITS_BData bf;     /* buffer of a temp block                     */
    int DiagOpt;  /* Option for diagonal inversion/solutiob     *
                   * opt =  1 -->> call luinv 
                   * opt == 2 -->> block inverted call dgemv    */

} ITS_VBILUSpar, *ITS_VBILUPtr; 

typedef struct ITS_ILUSpar_
{
    int n;
    ITS_CsPtr L;      /* L part elements                            */
    double *D;        /* diagonal elements                          */
    ITS_CsPtr U;      /* U part elements                            */
    int *work;        /* working buffer */

} ITS_ILUSpar, ITS_LDUmat, *ITS_ILUPtr;

typedef struct ITS_Per4Mat_ *ITS_P4Ptr;

/*------------------------------------------------------------
  | struct for storing the block LU factorization 
  | contains all the block factors except the 
  | data related to the last block. 
  | n       = size of current block
  | symperm = whether or not permutations are symmetric.
  |           used only in cleanP4..
  | nB      = size of B-block
  | L, U    = ILU factors of B-block
  | F, E    = sparse matrices in (1,2) and (2,1) 
  |           parts of matrix. 
  | perm    = (symmetric) permutation used in factorization
  |           comes from the independent set ordering
  | rperm   = unsymmetric permutation (rows) not used in this
  |           version -- but left here for compatibility..
  | D1, D2  = diagonal matrices (left, right) used for scaling
  |           if scaling option is turned on. Note that the 
  |           method works by scaling the whole matrix first
  |           (at any level) before anything else is done. 
  | wk     = a work vector of length n needed for various tasks
  |            [reduces number of calls to malloc]           
  |----------------------------------------------------------*/ 
typedef struct ITS_Per4Mat_
{
    int n;                  
    int nB; 
    int symperm;
    /*   LU factors  */
    ITS_SparMat *L;
    ITS_SparMat *U;
    /* E, F blocks   */
    ITS_SparMat *E;
    ITS_SparMat *F;
    int *rperm;       /* row permutation         */ 
    int *perm;        /* col. permutation        */ 
    double *D1 ;      /* diagonal scaling row    */  
    double *D2 ;      /* diagonal scaling columns*/  
    double *wk;       /* work array              */

    /* pointer to next and previous struct         */
    ITS_P4Ptr prev; 
    ITS_P4Ptr next;

} ITS_Per4Mat; 

/*------------------------------------------------------------
  | struct for storing data related to the last schur complement 
  | we need to store the C matrix associated with the last block
  | and the ILUT factorization of the related Schur complement.
  | 
  | n       = size of C block = size of Schur complement
  | C       = C block of last level matrix. 
  | L, U    = ILU factors of last schur complement. 
  |
  | meth[4] = parameters for defining variants in factorization 
  |           - see function readin for details
  | rperm    = row permutation used for very nonsymmetric matrices 
  |            [such as bottleneck transversal] -- NOT IN THIS VERSION
  | perm2     = unsymmetric permutation (columns) - used primarily
  |           for the ILUTP version of ILUT/.. 
  | D1, D2  = diagonal matrices (left, right) used for scaling
  |           if scaling option is turned on. Note that the 
  |           method works by scaling the whole matrix first
  |           (at any level) before anything else is done. 
  | wk     = a work vector of length n needed for various tasks
  |            [reduces number of calls to malloc]           
  |-----------------------------------------------------------*/
typedef struct ITS_ILUTFac *ITS_ILUTPtr;
typedef struct ITS_ILUTFac
{
    int n;                  

    /*-------------------- C matrix of last block */
    ITS_SparMat *C;

    /* LU factorization       */
    ITS_SparMat *L;
    ITS_SparMat *U;

    int *rperm;   /* row single-sinded permutation */
    int *perm;    /* column perm .                */
    int *perm2;   /* column permutation coming from pivoting in ILU */ 
    double *D1;
    double *D2;
    double *wk;

} ITS_ILUTSpar;

/* this is the arms preconditioner struct 
   | it consists of a linked list of p4mat structs
   | and the ILUt factorization (in the  form of an 
   | ITS_ILUTSpar struct  
   |---------------------------------------------- */
typedef struct ITS_ARMS_ *ITS_ARMSPtr;
typedef struct ITS_ARMS_
{
    int n;                   /* dimension of matrix */
    int nlev;                /* number of levels    */
    ITS_ILUTPtr ilus;        /* ILU for last level  */
    ITS_P4Ptr levmat;        /* level structure     */

} ITS_ARMSMat;

typedef struct ITS_CompressType
{
    int grp;   /* -1: begin new group, >=0: belong to grp-th row */
    int count; /* block size, valid only if grp = -1 */

} ITS_CompressType;

/*-------------------- 3 types of matrices so far */
typedef struct ITS_SMat
{
    int n; 
    int Mtype;           /*--  type 1 = CSR, 2 = VBCSR, 3 = LDU    */
    ITS_CsPtr CS;        /* place holder for a CSR/CSC type matrix */
    ITS_ILUPtr LDU;          /* struct for an LDU type matrix          */
    ITS_VBSPtr VBCSR;        /* place holder for a block matrix        */
    void (*matvec)(struct ITS_SMat*, double *, double *);

} ITS_SMat, *ITS_SMatptr;

/*-------------------- 3 types of matrices so far */
typedef struct ITS_SPre
{
    int Ptype;           /*-- Ptype 1 = ILU, 2 = VBILU, 3 = Crout */
    ITS_ILUPtr   ILU;        /* struct for an ILU type preconditioner */
    ITS_VBILUPtr VBILU;      /* struct for a block preconditioner */
    ITS_ARMSPtr ARMS;           /* struct for a block preconditioner */
    int (*precon) (double *, double *, struct ITS_SPre *); 

} ITS_SPre, *ITS_SPreptr;

typedef struct ITS_IOT_
{
    FILE *fout;                     /* output file handle              */
    char outfile[ITS_MAX_LINE];     /* output filename                 */
    char Fname[ITS_MAX_LINE];       /* full matrix path name           */
    char MatNam[ITS_MaxNamLen];     /* short name                      */
    char PrecMeth[ITS_MAX_LINE];    /* preconditioner being tested     */

    char type[4];               /* type for HB matrices            */
    int Fmt;                    /* matrix format type              */
    int ndim;                   /* matrix size                     */
    int nnz;                    /* number of nonzero               */

    /* parameters from inputs -----------------------------------------*/
    int im;                     /* Dim of Krylov subspace [fgmr]   */
    int maxits;                 /* maximum number of fgmres iters  */
    double tol;                 /* tolerance for stopping fgmres   */
    double eps;   /* for checking how close two rows of matrix are */
    int nparam;         /* number of tests for each preconditioner */
    int lfil0;                  /* initial lfil                    */
    int lfilInc;                /* increment for lfil              */
    double tol0;                /* initial drop tolerance          */
    double tolMul;              /* multiplier for tol              */
    int fill_lev;               /* initial level of fill for ILUK  */
    int fill_lev_inc;           /* increment for level of fill for ILUK */

    /* value always set to 1           */
    int perm_type;               /* indset perms (0) or PQ perms (1)*/

    /*                  or coarsen (2) */
    int Bsize;                   /* block size - dual role. see input file */

    /* for explanations */
    /* result for output ----------------------------------------------*/
    double rt_v;                /* compression rate of vertices    */
    double rt_e;                /* compression rate of edges       */
    double ceff;                /* compression efficiency          */
    double tm_h;                /* time for hash method  [vbilu]   */
    double tm_a;                /* time for angle method [vbilu]   */
    double tm_b;                /* time for initial blocks (s)     */
    double tm_p;                /* time for preconditioner (s)     */
    double tm_i;                /* time for iteration (s)          */
    double fillfact;            /* memory used during precondition */
    int its;                    /* number of iterations            */
    double enorm;               /* error norm:          || x- x0|| */
    double rnorm;               /* final residual norm: ||Ax-Ax0|| */

} ITS_IOT;

#endif  /* __VBLOCK_HEADER_H__ */
