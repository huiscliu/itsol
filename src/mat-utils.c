
#include "mat-utils.h"

#define TOL 1.e-17

/* *-------------------- inversion by svd
   This calls lapack routines for inverting a dense matrix.
   dgetrf and dgetri

   ON ENTRY
 ** A = square matrix of size n x n -- dimensioned with
 ** leading dimension =  n

 ON RETURN A contains the inverse  of the input matrix.
 */
int invGauss(int nn, double *A)
{
    int lWk, info;

    double *Wk;
    int *ipiv;

    lWk = 10 * nn;

    /*-------------------- trivial case nn = 1                     */
    if (nn == 1) {
        if (A[0] == 0.0)
            return 1;
        else {
            A[0] = 1.0 / A[0];
            return 0;
        }
    }
    /*-------------------- general case                              */

    Wk = (double *)malloc(lWk * sizeof(double));
    ipiv = (int *)malloc(nn * sizeof(int));
    if (Wk == NULL || ipiv == NULL)
        return -1;
    /*-------------------- get LU factorization with pivoting         */
    DGETRF(nn, nn, A, nn, ipiv, &info);
    if (info != 0)
        return info;
    /*-------------------- compute inverse                            */
    DGETRI(nn, A, nn, ipiv, Wk, lWk, &info);

    free(Wk);
    free(ipiv);
    return info;
}

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
int invSVD(int nn, double *A)
{
    int lWk, i, info;

    double *U, *VT, *S, *Wk;
    double tmp, nrm, one = 1.0, zero = 0.0;

    double tol = TOL;

    lWk = 5 * nn;

    U = (double *)malloc(nn * nn * sizeof(double));
    VT = (double *)malloc(nn * nn * sizeof(double));
    S = (double *)malloc(nn * sizeof(double));
    Wk = (double *)malloc(lWk * sizeof(double));

    if (U == NULL || VT == NULL || S == NULL || Wk == NULL)
        return -1;
    /*-------------------- trivial case nn = 1                     */
    if (nn == 1) {
        if (A[0] == 0.0)
            return 1;
        else {
            A[0] = one / A[0];
            return 0;
        }
    }
    /*-------------------- general case                              */
    dgesvd("A", "A", &nn, &nn, A, &nn, S, U, &nn, VT, &nn, Wk, &lWk, &info);
    if (S[0] == 0.0)
        return 1;
    nrm = S[0] * tol;
    /*-------------------- compute S\inv * VT                        */
    for (i = 0; i < nn; i++) {
        tmp = one / max(S[i], nrm);
        DSCAL(nn, tmp, &VT[i], nn);
    }
    /*-------------------- do [V^T S\inv ] * U^T                     */
    dgemm("t", "t", &nn, &nn, &nn, &one, VT, &nn, U, &nn, &zero, A, &nn);
    /*-------------------- Done -------------------------------------*/
    free(U);
    free(VT);
    free(S);
    free(Wk);
    return 0;
}
/*----------------------------------------------------------------------------
 * Diagonal scaling:
 * For the matrix with block diagonals D1, D2, ..., Dp :
 *       D1 x  x  x
 *       x  D2 x  x
 *  A =  x  x ... x
 *       x  x  x  Dp
 * simply take the block diagonal matrix
 *       D1 0  0  0
 *       0  D2 0  0
 *  D =  0  0 ... 0
 *       0  0  0  Dp
 * invert D and do A := inv(D)*A
 * then the diagonal blocks of A are now identities.
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * vbmat    = block matrix stored in VBSpaFmt format -- see type-defs.h for
 *            details on format, the block sizes might be different
 * on return:
 * ==========
 * vbmat    = inv(D)*vbmat
 *--------------------------------------------------------------------------*/
int diag_scal(vbsptr vbmat)
{
    int i, j, k, dim, sz, size, ierr = 0, col;
    double one = 1.0, zero = 0.0;
    int nzcount, n = vbmat->n, *bsz = vbmat->bsz, *ja;
    int bufsz = sizeof(double) * MAX_BLOCK_SIZE * MAX_BLOCK_SIZE;
    BData *ba, *D, buf;

    D = (BData *) Malloc(sizeof(BData) * n, "diag_scal");
    buf = (BData) Malloc(bufsz, "diag_scal");

    for (i = 0; i < n; i++) {
        nzcount = vbmat->nzcount[i];
        ja = vbmat->ja[i];

        for (j = 0; j < nzcount; j++) {
            if (ja[j] != i)
                continue;

            dim = B_DIM(bsz, i);
            size = sizeof(double) * dim * dim;
            D[i] = (BData) Malloc(size, "diag_scal");
            memcpy(D[i], vbmat->ba[i][j], size);

            ierr = invSVD(dim, D[i]);
            if (ierr != 0) {
                for (k = 0; k < i; k++)
                    free(D[k]);
                free(D);
                fprintf(stderr, "error: Singular diagonal block...\n");
                return -2;
            }
        }
    }

    for (i = 0; i < n; i++) {
        dim = B_DIM(bsz, i);
        nzcount = vbmat->nzcount[i];
        ja = vbmat->ja[i];
        ba = vbmat->ba[i];

        for (j = 0; j < nzcount; j++) {
            col = ja[j];
            sz = B_DIM(bsz, col);
            DGEMM("n", "n", dim, sz, dim, one, D[i], dim, ba[j], dim, zero, buf, dim);
            copyBData(dim, sz, ba[j], buf, 0);
        }
    }

    vbmat->D = D;
    free(buf);
    return 0;
}

/*---------------------------------------------------------------------
  | This function does y = inv(D) x, where D is diagonals of A.
  |----------------------------------------------------------------------
  | on entry:
  | vbmat = the matrix (in BSpaFmt form)
  | x     = a vector
  |
  | on return
  | y     = the product inv(D) * x
  |--------------------------------------------------------------------*/
int diagvec(vbsptr vbmat, BData x, BData y)
{
    int i, n = vbmat->n, *bsz = vbmat->bsz, dim, sz = 1;
    double zero = 0.0, one = 1.0;
    BData *D = vbmat->D;

    for (i = 0; i < n; i++) {
        dim = B_DIM(bsz, i);
        DGEMM("n", "n", dim, sz, dim, one, D[i], dim, x + bsz[i], dim, zero, y + bsz[i], dim);
    }
    return 0;
}

/*---------------------------------------------------------------------
  | This function does the matrix vector product y = A x.
  |----------------------------------------------------------------------
  | on entry:
  | mata  = the matrix (in SpaFmt form)
  | x     = a vector
  |
  | on return
  | y     = the product A * x
  |--------------------------------------------------------------------*/
void matvec(csptr mata, double *x, double *y)
{
    int i, k, *ki;
    double *kr;

    for (i = 0; i < mata->n; i++) {
        y[i] = 0.0;
        kr = mata->ma[i];
        ki = mata->ja[i];

        for (k = 0; k < mata->nzcount[i]; k++)
            y[i] += kr[k] * x[ki[k]];
    }
}

void vbmatvec(vbsptr vbmat, double *x, double *y)
{
    int i, j, nzcount, col, inc = 1, dim, sz, nBs, nBsj;
    int n = vbmat->n, *ja, *bsz = vbmat->bsz;
    double one = 1.0;
    BData *ba;

    for (i = 0; i < n; i++) {
        nBs = bsz[i];
        dim = B_DIM(bsz, i);

        for (j = 0; j < dim; j++)
            y[nBs + j] = 0;

        nzcount = vbmat->nzcount[i];
        ja = vbmat->ja[i];
        ba = vbmat->ba[i];
        for (j = 0; j < nzcount; j++) {
            col = ja[j];
            nBsj = bsz[col];
            sz = B_DIM(bsz, col);

            DGEMV("n", dim, sz, one, ba[j], dim, &x[nBsj], inc, one, &y[nBs], inc);
        }
    }
}

/*---------------------------------------------------------------------
  | This function does the forward solve L x = b.
  | Can be done in place.
  |----------------------------------------------------------------------
  | on entry:
  | mata  = the matrix (in SpaFmt form)
  | b     = a vector
  |
  | on return
  | x     = the solution of L x = b 
  |--------------------------------------------------------------------*/
void Lsol(csptr mata, double *b, double *x)
{
    int i, k;
    double *kr;
    int *ki;

    for (i = 0; i < mata->n; i++) {
        x[i] = b[i];

        if (mata->nzcount[i] > 0) {
            kr = mata->ma[i];
            ki = mata->ja[i];

            for (k = 0; k < mata->nzcount[i]; k++)
                x[i] -= kr[k] * x[ki[k]];
        }
    }
}

/*---------------------------------------------------------------------
  | This function does the backward solve U x = b.
  | Can be done in place.
  |----------------------------------------------------------------------
  | on entry:
  | mata  = the matrix (in SpaFmt form)
  | b    = a vector
  |
  | on return
  | x     = the solution of U * x = b 
  |
  |---------------------------------------------------------------------*/
void Usol(csptr mata, double *b, double *x)
{
    int i, k, *ki;
    double *kr;

    for (i = mata->n - 1; i >= 0; i--) {
        kr = mata->ma[i];
        ki = mata->ja[i];
        x[i] = b[i];

        for (k = 1; k < mata->nzcount[i]; k++)
            x[i] -= kr[k] * x[ki[k]];

        x[i] *= kr[0];
    }
}

/*---------------------------------------------------------------------
  | This function does the (block) forward elimination in ARMS
  |                       new       old
  |     |            |  |     |    |    |
  |     | L        0 |  | wx1 |    | x1 |
  |     |            |  |     | =  |    | 
  |     | EU^{-1}  I |  | wx2 |    | x2 |
  |     |            |  |     |    |    |
  | x used and not touched -- or can be the same as wk.
  |--------------------------------------------------------------------*/
int descend(p4ptr levmat, double *x, double *wk)
{
    /*  local variables   */
    int j, len = levmat->n, lenB = levmat->nB, *iperm = levmat->rperm;
    double *work = levmat->wk;

    for (j = 0; j < len; j++)
        work[iperm[j]] = x[j];

    Lsol(levmat->L, work, wk);  /* sol:   L x = x                 */
    Usol(levmat->U, wk, work);  /* sol:   U work(2) = work         */

    /*-------------------- compute x[lenb:.] = x [lenb:.] - E * work(1) */
    matvecz(levmat->E, work, &work[lenB], &wk[lenB]);
    return 0;
}

/*---------------------------------------------------------------------
  | This function does the (block) backward substitution: 
  |
  |     |            |  |     |    |    |
  |     | U  L^{-1}F |  | wk1 |    | x1 |
  |     |            |  |     | =  |    |
  |     | 0       S  |  | wk2 |    | x2 |  <<-- x2 already computed.
  |     |            |  |     |    |    |       and we need x1
  |
  |    with x2 = S^{-1} wk2 [assumed to have been computed ] 
  |--------------------------------------------------------------------*/
int ascend(p4ptr levmat, double *x, double *wk)
{
    int j, len = levmat->n, lenB = levmat->nB, *qperm = levmat->perm;
    double *work = levmat->wk;

    matvec(levmat->F, &x[lenB], work);  /*  work = F * x_2   */
    Lsol(levmat->L, work, work);        /*  work = L \ work    */

    for (j = 0; j < lenB; j++)  /*  wk1 = wk1 - work  */
        work[j] = x[j] - work[j];

    Usol(levmat->U, work, work);        /*  wk1 = U \ wk1 */
    memcpy(&work[lenB], &x[lenB], (len - lenB) * sizeof(double));

    /*---------------------------------------
      |   apply reverse permutation
      |--------------------------------------*/
    for (j = 0; j < len; j++)
        wk[j] = work[qperm[j]];

    return 0;
}

/*---------------------------------------------------------------------
  | This function does the matrix vector  z = y - A x.
  |----------------------------------------------------------------------
  | on entry:
  | mata  = the matrix (in SpaFmt form)
  | x, y   = two input vector
  |
  | on return
  | z    = the result:  y - A * x
  | z-location must be different from that of x 
  | i.e., y and x are used but not modified.
  |--------------------------------------------------------------------*/
void matvecz(csptr mata, double *x, double *y, double *z)
{
    int i, k, *ki;
    double *kr, t;

    for (i = 0; i < mata->n; i++) {
        kr = mata->ma[i];
        ki = mata->ja[i];
        t = y[i];

        for (k = 0; k < mata->nzcount[i]; k++)
            t -= kr[k] * x[ki[k]];

        z[i] = t;
    }
}

/* Macro L-solve -- corresponds to left (L) part of arms
   |  preconditioning operation -- 
   |  on entry : 
   |   x =  right- hand side to be operated on by the preconditioner
   |  on return : x is overwritten
   |   x =  output result of operation 
   |  
   |  Note : in-place operation -- b and x can occupy the same space..
   | --------------------------------------------------------------------*/
p4ptr Lvsol2(double *x, int nlev, p4ptr levmat, ilutptr ilusch)
{
    int nloc = levmat->n, first, lenB;
    p4ptr last = levmat;

    /*-------------------- take care of  special cases :  nlev==0 --> lusol  */
    if (nlev == 0) {
        SchLsol(ilusch, x);
        return (last);
    }

    first = 0;
    /*-------------------- descend                                      */
    while (levmat) {
        nloc = levmat->n;
        lenB = levmat->nB;
        /*-------------------- left scaling                                  */
        if (levmat->D1 != NULL)
            dscale(nloc, levmat->D1, &x[first], &x[first]);
        /*--------------------  RESTRICTION/ DESCENT OPERATION  */
        if (lenB)
            descend(levmat, &x[first], &x[first]);
        first += lenB;
        last = levmat;
        levmat = levmat->next;
    }

    SchLsol(ilusch, &x[first]);

    return last;
}

/* Macro U-solve -- corresponds to right (U) part of arms
   |  preconditioning operation -- 
   |  on entry : 
   |  b  =  right- hand side to be operated on by the preconditioner
   |  on return  = x has been overwritten =
   |  x  =  output result of operation 
   |  
   |  Note : in-place operation -- b and x  can occupy the same space..
   | --------------------------------------------------------------------*/
int Uvsol2(double *x, int nlev, int n, p4ptr levmat, ilutptr ilusch)
{
    int nloc, lenB, first;
    if (nlev == 0) {
        SchUsol(ilusch, x);
        return (0);
    }

    /*-------------------- general case                               */
    nloc = levmat->n;
    lenB = levmat->nB;
    first = n - nloc;
    /*-------------------- last level                                 */
    first += lenB;
    SchUsol(ilusch, &x[first]);
    /*-------------------- other levels                               */
    while (levmat) {
        nloc = levmat->n;
        first -= levmat->nB;
        if (levmat->n)
            ascend(levmat, &x[first], &x[first]);
        /*-------------------- right scaling */
        if (levmat->D2 != NULL)
            dscale(nloc, levmat->D2, &x[first], &x[first]);
        levmat = levmat->prev;
    }
    return 0;
}

/* combined preconditioning operation -- combines the
   |  left and right actions. 
   | 
   |  on entry : 
   |   x =  right- hand side to be operated on by the preconditioner
   |  on return : x is overwritten - 
   |   x =  output result of operation 
   |  
   |  Note : in-place operation -- b and x can occupy the same space..
   | --------------------------------------------------------------------*/
int armsol2(double *x, arms Prec)
{
    p4ptr levmat = Prec->levmat;
    ilutptr ilusch = Prec->ilus;
    int nlev = Prec->nlev;
    int n = levmat->n;
    p4ptr last;

    if (nlev == 0) {
        n = ilusch->n;
        SchLsol(ilusch, x);
        SchUsol(ilusch, x);
        return 0;
    }

    last = Lvsol2(x, nlev, levmat, ilusch);
    Uvsol2(x, nlev, n, last, ilusch);

    return 0;
}

/*---------------------------------------------------------------------
  |  Forward solve for Schur complement part = 
  |----------------------------------------------------------------------
  | on entry:
  | ilusch  = the LU matrix as provided from the ILU functions.
  | y       = the right-hand-side vector
  |
  | on return
  | y       = solution of LU x = y. [overwritten] 
  |---------------------------------------------------------------------*/
void SchLsol(ilutptr ilusch, double *y)
{
    int n = ilusch->n, j, *perm = ilusch->rperm;
    double *work = ilusch->wk;

    /*-------------------- begin: right scaling                          */
    if (ilusch->D1 != NULL)
        dscale(n, ilusch->D1, y, y);
    /*-------------------- ONE SIDED ROW PERMS */
    if (perm != NULL) {
        for (j = 0; j < n; j++)
            work[perm[j]] = y[j];
        /*--------------------  L solve proper */
        Lsol(ilusch->L, work, y);
    }
    else
        Lsol(ilusch->L, y, y);
}

/*---------------------------------------------------------------------
  | U-solve for Schur complement  - 
  |----------------------------------------------------------------------
  | on entry:
  | ilusch  = the LU matrix as provided from the ILU functions.
  | y       = the right-hand-side vector
  |
  | on return 
  | y       = solution of U x = y. [overwritten on y] 
  |----------------------------------------------------------------------*/
void SchUsol(ilutptr ilusch, double *y)
{
    int n = ilusch->n, j, *perm = ilusch->perm, *cperm;
    double *work = ilusch->wk;
    /* -------------------- begin by U-solving */
    /*-------------------- CASE: column pivoting  used (as in ILUTP) */
    if (ilusch->perm2 != NULL) {
        Usol(ilusch->U, y, y);
        cperm = ilusch->perm2;
        for (j = 0; j < n; j++)
            work[cperm[j]] = y[j];
    }
    else
        /*-------------------- CASE: no column pivoting  used                   */
        Usol(ilusch->U, y, work);
    /*-------------------- generic permutation                              */
    if (perm != NULL) {
        for (j = 0; j < n; j++)
            y[j] = work[perm[j]];
    }
    else
        memcpy(y, work, n * sizeof(double));

    /*-------------------- case when diagonal scaling is done on columns    */
    if (ilusch->D2 != NULL)
        dscale(n, ilusch->D2, y, y);
}

/*--------------------------------------------------------
 *    does the operation y = inv(a) * x
 *    where a has already been factored by Gauss.
 *    LUy = x
 *------------------------------------------------------*/
void luinv(int n, double *a, double *x, double *y)
{
    int i, j, bsA, bsB;
    double sum;

    /* Ly0 = x -- use Lsol ? */
    for (i = 0; i < n; i++) {
        sum = x[i];
        bsA = i - n;
        for (j = 0; j < i; j++) {
            bsA += n;
            sum -= a[bsA] * y[j];       /* a(i,j) * y(j) */
        }
        y[i] = sum;
    }

    /* Uy = y0 */
    bsB = i * n;
    for (i = n - 1; i >= 0; i--) {
        sum = y[i];
        bsB -= n;
        bsA = i + bsB;
        for (j = i + 1; j < n; j++) {
            bsA += n;
            sum -= a[bsA] * y[j];       /* a(i,j) * y(j) */
        }
        y[i] = sum * a[bsB + i];        /* a(i,i) */
    }
}

/*----------------------------------------------------------------------
 *    performs a forward followed by a backward solve
 *    for LU matrix as produced by iluk
 *    y  = right-hand-side 
 *    x  = solution on return 
 *    lu = LU matrix as produced by iluk. 
 *--------------------------------------------------------------------*/
int lusolC(double *y, double *x, iluptr lu)
{
    int n = lu->n, i, j, nzcount, *ja;
    double *D;
    csptr L, U;

    L = lu->L;
    U = lu->U;
    D = lu->D;

    /* Block L solve */
    for (i = 0; i < n; i++) {
        x[i] = y[i];
        nzcount = L->nzcount[i];
        ja = L->ja[i];
        for (j = 0; j < nzcount; j++) {
            x[i] -= x[ja[j]] * L->ma[i][j];
        }
    }
    /* Block -- U solve */
    for (i = n - 1; i >= 0; i--) {
        nzcount = U->nzcount[i];
        ja = U->ja[i];
        for (j = 0; j < nzcount; j++) {
            x[i] -= x[ja[j]] * U->ma[i][j];
        }
        x[i] *= D[i];
    }
    return (0);
}

/*----------------------------------------------------------------------
 *    performs a forward followed by a backward solve
 *    for LU matrix as produced by iluc
 *    y  = right-hand-side
 *    x  = solution on return
 *    lu = LU matrix as produced by iluc.
 *--------------------------------------------------------------------*/
int lumsolC(double *y, double *x, iluptr lu)
{
    int n = lu->n, i, j, nzcount, nnzL, *ia, *ja;
    double *D = lu->D, *ma;
    csptr L = lu->L;
    csptr U = lu->U;

    for (i = 0; i < n; i++)
        x[i] = y[i];
    /*-------------------- L solve */
    for (i = 0; i < n; i++) {
        nnzL = L->nzcount[i];
        ia = L->ja[i];
        ma = L->ma[i];
        for (j = 0; j < nnzL; j++) {
            x[ia[j]] -= ma[j] * x[i];
        }
    }
    /*-------------------- U solve */
    for (i = n - 1; i >= 0; i--) {
        nzcount = U->nzcount[i];
        ja = U->ja[i];
        ma = U->ma[i];
        for (j = 0; j < nzcount; j++) {
            x[i] -= ma[j] * x[ja[j]];
        }
        x[i] *= D[i];
    }

    return 0;
}

/*----------------------------------------------------------------------
 *    performs a forward followed by a backward block solve
 *    for LU matrix as produced by VBILUT
 *    y  = right-hand-side 
 *    x  = solution on return 
 *    lu = LU matrix as produced by VBILUT
 *
 *    note: lu->bf is used to store vector
 *--------------------------------------------------------------------*/
int vblusolC(double *y, double *x, vbiluptr lu)
{
    int n = lu->n, *bsz = lu->bsz, i, j, bi, icol, dim, sz;
    int nzcount, nBs, nID, *ja, inc = 1, OPT;
    double *data, alpha = -1.0, beta = 1.0, alpha2 = 1.0, beta2 = 0.0;
    vbsptr L, U;
    BData *D, *ba;

    L = lu->L;
    U = lu->U;
    D = lu->D;
    OPT = lu->DiagOpt;
    /* Block L solve */
    for (i = 0; i < n; i++) {
        dim = B_DIM(bsz, i);
        nBs = bsz[i];
        for (j = 0; j < dim; j++) {
            nID = nBs + j;
            x[nID] = y[nID];
        }

        nzcount = L->nzcount[i];
        ja = L->ja[i];
        ba = L->ba[i];
        for (j = 0; j < nzcount; j++) {
            icol = ja[j];
            sz = B_DIM(bsz, icol);
            data = ba[j];
            DGEMV("n", dim, sz, alpha, data, dim, x + bsz[icol], inc, beta, x + nBs, inc);
        }
    }
    /* Block -- U solve */
    for (i = n - 1; i >= 0; i--) {
        dim = B_DIM(bsz, i);
        nzcount = U->nzcount[i];
        nBs = bsz[i];
        ja = U->ja[i];
        ba = U->ba[i];
        for (j = 0; j < nzcount; j++) {
            icol = ja[j];
            sz = B_DIM(bsz, icol);
            data = ba[j];
            DGEMV("n", dim, sz, alpha, data, dim, x + bsz[icol], inc, beta, x + nBs, inc);
        }
        data = D[i];
        if (OPT == 1)
            luinv(dim, data, x + nBs, lu->bf);
        else
            DGEMV("n", dim, dim, alpha2, data, dim, x + nBs, inc, beta2, lu->bf, inc);

        for (bi = 0; bi < dim; bi++) {
            x[nBs + bi] = lu->bf[bi];
        }
    }

    return 0;
}

/*----------------------------------------------------------------------
  |
  | This subroutine permutes the rows of a matrix in SpaFmt format. 
  | rperm  computes B = P A  where P is a permutation matrix.  
  | The permutation P is defined through the array perm: for each j, 
  | perm[j] represents the destination row number of row number j. 
  |
  |-----------------------------------------------------------------------
  | on entry:
  |----------
  | (amat) = a matrix stored in SpaFmt format.
  |
  |
  | on return:
  | ----------
  | (amat) = P A stored in SpaFmt format.
  |
  | integer value returned:
  |             0   --> successful return.
  |             1   --> memory allocation error.
  |---------------------------------------------------------------------*/
int rpermC(csptr mat, int *perm)
{
    int **addj, *nnz, i, size = mat->n;
    double **addm;
    addj = (int **)Malloc(size * sizeof(int *), "rpermC");
    addm = (double **)Malloc(size * sizeof(double *), "rpermC");
    nnz = (int *)Malloc(size * sizeof(int), "rpermC");
    for (i = 0; i < size; i++) {
        addj[perm[i]] = mat->ja[i];
        addm[perm[i]] = mat->ma[i];
        nnz[perm[i]] = mat->nzcount[i];
    }
    for (i = 0; i < size; i++) {
        mat->ja[i] = addj[i];
        mat->ma[i] = addm[i];
        mat->nzcount[i] = nnz[i];
    }
    free(addj);
    free(addm);
    free(nnz);
    return 0;
}

/*----------------------------------------------------------------------
  |
  | This subroutine permutes the columns of a matrix in SpaFmt format.
  | cperm computes B = A P, where P is a permutation matrix.
  | that maps column j into column perm(j), i.e., on return 
  | The permutation P is defined through the array perm: for each j, 
  | perm[j] represents the destination column number of column number j. 
  |
  |-----------------------------------------------------------------------
  | on entry:
  |----------
  | (mat) = a matrix stored in SpaFmt format.
  |
  |
  | on return:
  | ----------
  | (mat) = A P stored in SpaFmt format.
  |
  | integer value returned:
  |             0   --> successful return.
  |             1   --> memory allocation error.
  |---------------------------------------------------------------------*/
int cpermC(csptr mat, int *perm)
{
    int i, j, *newj, size = mat->n, *aja;

    newj = (int *)Malloc(size * sizeof(int), "cpermC");
    for (i = 0; i < size; i++) {
        aja = mat->ja[i];
        for (j = 0; j < mat->nzcount[i]; j++)
            newj[j] = perm[aja[j]];

        for (j = 0; j < mat->nzcount[i]; j++)
            aja[j] = newj[j];
        mat->ja[i] = aja;
    }
    free(newj);
    return 0;
}

/*----------------------------------------------------------------------
  |
  | This subroutine permutes the rows and columns of a matrix in 
  | SpaFmt format.  dperm computes B = P^T A P, where P is a permutation 
  | matrix.
  |
  |-----------------------------------------------------------------------
  | on entry:
  |----------
  | (amat) = a matrix stored in SpaFmt format.
  |
  |
  | on return:
  | ----------
  | (amat) = P^T A P stored in SpaFmt format.
  |
  | integer value returned:
  |             0   --> successful return.
  |             1   --> memory allocation error.
  |---------------------------------------------------------------------*/
int dpermC(csptr mat, int *perm)
{
    if (rpermC(mat, perm))
        return 1;
    if (cpermC(mat, perm))
        return 1;
    return 0;
}

/*----------------------------------------------------------------------
  | Finds the compressed transpose of a matrix stored in SpaFmt format.
  | Patterns only.
  |-----------------------------------------------------------------------
  | on entry:
  |----------
  | (amat)     = a matrix stored in SpaFmt format.
  | (compress) = quotient graph of matrix amat
  |
  | on return:
  | ----------
  | (bmat)     = the compressed transpose of (mata) stored in SpaFmt
  |              format.
  |
  | integer value returned:
  |             0   --> successful return.
  |             1   --> memory allocation error.
  |---------------------------------------------------------------------*/
int CSparTran(csptr amat, csptr bmat, CompressType * compress)
{
    int i, j, *ind, nzcount, pos, size = amat->n, *aja;
    ind = bmat->nzcount;

    for (i = 0; i < size; i++)
        ind[i] = 0;
    /*-------------------- compute lengths  */
    for (i = 0; i < size; i++) {
        if (compress[i].grp != -1)
            continue;
        aja = amat->ja[i];
        nzcount = amat->nzcount[i];
        for (j = 0; j < nzcount; j++) {
            pos = aja[j];
            if (compress[pos].grp == -1) {
                ind[pos]++;
            }
        }
    }

    /*--------------------  allocate space  */
    for (i = 0; i < size; i++) {
        if (ind[i] == 0) {
            bmat->ja[i] = NULL;
            continue;
        }
        bmat->ja[i] = (int *)Malloc(ind[i] * sizeof(int), "CSparTran");
        ind[i] = 0;             /* indicate next available position of each row */
    }
    /*--------------------  now do the actual copying  */
    for (i = 0; i < size; i++) {
        if (compress[i].grp != -1)
            continue;
        aja = amat->ja[i];
        nzcount = amat->nzcount[i];
        for (j = 0; j < nzcount; j++) {
            pos = aja[j];
            if (compress[pos].grp == -1) {
                bmat->ja[pos][ind[pos]] = i;
                ind[pos]++;
            }
        }
    }
    return 0;
}

double vbnorm2(int sz, double *a)
{
    int tmp = 1;

    return DNRM2(sz, a, tmp) / (double)sz;
}

/*-----------------------------------------------------------*/
int condestLU(iluptr lu, FILE * fp)
{
    int n = lu->n, i;
    double norm = 0.0;
    double *y = (double *)Malloc(n * sizeof(double), "condestLU");
    double *x = (double *)Malloc(n * sizeof(double), "condestLU");

    for (i = 0; i < n; i++)
        y[i] = 1.0;
    lusolC(y, x, lu);
    for (i = 0; i < n; i++)
        norm = max(norm, fabs(x[i]));
    fprintf(fp, "ILU inf-norm lower bound : %16.2f\n", norm);
    free(y);
    free(x);
    if (norm > 1e30)
        return -1;
    return 0;
}

/*-------------------- simple estimate of cond. number of precon */
int condestArms(arms armspre, double *y, FILE * fp)
{
    int n = armspre->n, i;
    double norm = 0.0;

    for (i = 0; i < n; i++)
        y[i] = 1.0;
    armsol2(y, armspre);
    for (i = 0; i < n; i++) {
        norm = max(norm, fabs(y[i]));
    }
    fprintf(fp, "ARMS inf-norm lower bound = : %16.2f\n", norm);
    if (norm > 1e30) {
        return -1;
    }
    return 0;
}

int VBcondestC(vbiluptr lu, FILE * fp)
{
    int n = lu->n, i, ndim = lu->bsz[n];
    double norm = 0.0;
    double *y = (double *)Malloc(ndim * sizeof(double), "condestLU");
    double *x = (double *)Malloc(ndim * sizeof(double), "condestLU");

    for (i = 0; i < ndim; i++)
        y[i] = 1.0;
    vblusolC(y, x, lu);
    for (i = 0; i < ndim; i++)
        norm = max(norm, fabs(x[i]));

    fprintf(fp, "VBILU inf-norm lower bound : %16.2f\n", norm);
    free(y);
    free(x);
    if (norm > 1e30) {
        return -1;
    }
    return 0;
}

/*-------------------------------------------------------------------
  | This function does the matrix vector product y = A x. COL storage
  |--------------------------------------------------------------------
  | on entry:
  | mat  = the matrix (in SpaFmt form -- COLUMN) | x = a vector
  | on return
  | y     = the product A * x
  |--------------------------------------------------------------------*/
void matvecC(csptr mat, double *x, double *y)
{
    int n = mat->n, i, k, *ki;
    double *kr;

    for (i = 0; i < n; i++)
        y[i] = 0.0;

    for (i = 0; i < n; i++) {
        kr = mat->ma[i];
        ki = mat->ja[i];

        for (k = 0; k < mat->nzcount[i]; k++)
            y[ki[k]] += kr[k] * x[i];
    }
}

void matvecCSR(SMatptr mat, double *x, double *y)
{
    matvec(mat->CS, x, y);
}

void matvecCSC(SMatptr mat, double *x, double *y)
{
    matvecC(mat->CS, x, y);
}

void matvecVBR(SMatptr mat, double *x, double *y)
{
    vbmatvec(mat->VBCSR, x, y);
}

int preconILU(double *x, double *y, SPreptr mat)
{
    /*-------------------- precon for csr format using the SPre struct*/
    return lusolC(x, y, mat->ILU);
}

int preconVBR(double *x, double *y, SPreptr mat)
{
    /*-------------------- precon for ldu format using the SPre struct*/
    return vblusolC(x, y, mat->VBILU);
}

int preconLDU(double *x, double *y, SPreptr mat)
{
    /*-------------------- precon for vbr format using the SPre struct*/
    return lumsolC(x, y, mat->ILU);
}

int preconARMS(double *x, double *y, SPreptr mat)
{
    /*-------------------- precon for ldu format using the SPre struct*/
    int n = (mat->ARMS)->n;
    memcpy(y, x, n * sizeof(double));
    return armsol2(y, mat->ARMS);
}

typedef struct __KeyType {
    int var;                    /* row number */
    int key;                    /* hash value */

} KeyType;

static int KeyComp(const void *vfst, const void *vsnd)
{
    KeyType *fst = (KeyType *) vfst, *snd = (KeyType *) vsnd;
    if (fst->key == snd->key) {
        if (fst->var < snd->var)
            return -1;
        return 1;
    }
    if (fst->key < snd->key)
        return -1;
    return 1;
}

/*----------------------------------------------------------------------------
 * Setup Blocks ( rows and columns might be permuted to get better results )
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
 *        if cos( <row_i, row_j> ) = (row_i,row_j)/|row_i||row_j| is > eps,
 *        we merge row_i and row_j by resetting
 *        group[j] = i and size[i] = size[i]+size[j]
 *--------------------------------------------------------------------------*/
int init_blocks(csptr csmat, int *pnBlock, int **pnB, int **pperm, double eps, double *t_hash, double *t_angle)
{
    int n = csmat->n, nBlock = 0, i, j, k;
    csptr at = NULL;
    KeyType *group = NULL;
    CompressType *compress = NULL;
    int *perm = NULL, *nB = NULL;
    int nzcount0, nzcount, key0, key, *ja0, *ja, row0, row, newblock;
    int *iw = NULL, *jbuf = NULL;
    int cnt, pos, nnz_i, row_j, col, bkcnt;
    int nextBlockID, nextBlockPos, belongTo, grp;
    double eps_2 = eps * eps, t1, t2;

    t1 = sys_timer();           /* begin Hash method timer */
    group = (KeyType *) Malloc(n * sizeof(KeyType), "init_blocks");
    compress = (CompressType *) Malloc(n * sizeof(CompressType), "init_blocks");
    perm = (int *)Malloc(n * sizeof(int), "init_blocks");
    iw = perm;                  /* iw and perm array can share memory here because they will
                                 * never be used at the same time */
    for (i = 0; i < n; i++) {
        iw[i] = 0;
        compress[i].grp = -1;
    }
    /*-------------------- compress matrix based on hash algorithm */
    /*-------------------- get hash value of each row */
    for (i = 0; i < n; i++) {
        nzcount = csmat->nzcount[i];
        key = 0;
        ja = csmat->ja[i];
        for (j = 0; j < nzcount; j++)
            key += ja[j] + 1;
        group[i].key = key;
        group[i].var = i;
    }
    /*-------------------- sort rows -- uses function KeyComp */
    qsort(group, n, sizeof(KeyType), KeyComp);

    /*-------------------- compress matrix */
    for (i = 0; i < n; i++) {
        row0 = group[i].var;
        if (compress[row0].grp != -1)
            continue;           /* already assigned */
        key0 = group[i].key;
        nzcount0 = csmat->nzcount[row0];
        ja0 = csmat->ja[row0];
        /*-------------------- beginning of new block. set .grp and .count */
        compress[row0].grp = -1;
        compress[row0].count = 1;
        /*-------------------- loop over all rows having same check-sum keys */
        for (j = i + 1; j < n; j++) {
            key = group[j].key;
            if (key != key0)
                break;
            row = group[j].var;
            if (compress[row].grp != -1)
                continue;       /* already assigned */
            nzcount = csmat->nzcount[row];
            if (nzcount != nzcount0)
                continue;
            ja = csmat->ja[row];
            newblock = 0;
            /*-------------------- compare patterns of the rows             */
            for (k = 0; k < nzcount; k++)
                iw[ja0[k]] = 1;
            for (k = 0; k < nzcount; k++) {
                if (iw[ja[k]] == 0) {
                    newblock = 1;
                    break;
                }
            }
            for (k = 0; k < nzcount; k++)
                iw[ja0[k]] = 0; /* reset iw */
            /*-------------------- row belongs to group row0                    */
            if (!newblock) {
                compress[row].grp = row0;
                compress[row0].count++;
            }
        }
    }
    t2 = sys_timer();           /* end Hash method timer */
    *t_hash = t2 - t1;

    t1 = sys_timer();           /* begin angle method timer */
    nB = (int *)Malloc(n * sizeof(int), "init_blocks");
    jbuf = (int *)Malloc(n * sizeof(int), "init_blocks");

    /*-------------------- compress matrix based on angle algorithm */
    /*-------------------- calculate compressed A^T                 */
    at = (csptr) Malloc(sizeof(SparMat), "init_blocks");
    setupCS(at, n, 0);
    if (CSparTran(csmat, at, compress) != 0)
        return -1;

    /*----------------------------------------------------------------------------
     * only the row representing beginning of block satisfies:
     *    compress[row].grp = -1, so far.
     * how many such rows is up to the compression rate of hash compression
     * algorithm we did above
     *--------------------------------------------------------------------------*/

    /*---------------------------------------------------------------
     * use group to backup original compressed matrix by Hash method.
     * It is very important because the array 'compress' will be changed
     * during Angle method. Or, we'll get incorrect inner product.
     *--------------------------------------------------------------*/
    for (i = 0; i < n; i++) {
        group[i].var = compress[i].grp;
        group[i].key = compress[i].count;
    }

    for (i = 0; i < n; i++) {
        if (compress[i].grp != -1)
            continue;
        nB[nBlock] = compress[i].count; /* !!! not 1 here */
        cnt = 0;
        /*-------------------- calculate (u,v_j ), j = i+1,...,n, using product
         *-------------------- algorithm of A * A_T */
        nnz_i = csmat->nzcount[i];
        for (j = 0; j < nnz_i; j++) {
            row_j = csmat->ja[i][j];
            if (group[row_j].var != -1) /* i.e. original compress[row_j].grp */
                continue;
            bkcnt = group[row_j].key;   /* i.e. original compress[row_j].count */
            for (k = at->nzcount[row_j] - 1; k >= 0; k--) {
                col = at->ja[row_j][k];
                if (col <= i)
                    break;
                if (compress[col].grp != -1)
                    continue;   /* needed because compress
                                   array is dynamically updated */
                if (iw[col] == 0) {     /* new nonzero of (u,v_j) */
                    jbuf[cnt] = col;
                    cnt++;
                }
                iw[col] += bkcnt;       /* correct for matrix with symmetric pattern */
            }
        }
        /*-------------------- set group for row i and reset iw */
        for (j = 0; j < cnt; j++) {
            pos = jbuf[j];
            if (iw[pos] * iw[pos] >= eps_2 * nnz_i * csmat->nzcount[pos]) {
                compress[pos].grp = i;
                nB[nBlock] += compress[pos].count;      /* !!! not 1 here */
            }
            iw[pos] = 0;        /* reset iw */
        }
        nBlock++;               /* begin new block, add block count by 1 */
    }                           /* end loop i */

    /*-------------------- free group                                   */
    if (group) {
        /* no need group array any more */
        free(group);
        group = NULL;
    }

    *pnBlock = nBlock;
    *pnB = (int *)Malloc(nBlock * sizeof(int), "init_blocks");
    for (i = 0; i < nBlock; i++) {
        if (nB[i] > MAX_BLOCK_SIZE) {
            fprintf(stderr, "Block of size = %d exceeds MAX_BLOCK_SIZE\n", nB[i]);
            return -1;
        }
        (*pnB)[i] = nB[i];
    }

    /*-------------------- calculate permutation array -  Array nB will
     * be used to store next  available position in each  block */
    nextBlockID = 0;
    nextBlockPos = 0;
    for (i = 0; i < n; i++) {
        if (compress[i].grp == -1) {
            perm[i] = nextBlockPos;
            nextBlockPos += (*pnB)[nextBlockID++];
            nB[i] = 1;
        }
        else {
            belongTo = compress[i].grp;
            grp = compress[belongTo].grp;
            if (grp != -1)      /* should find the final beginning of block */
                belongTo = grp;
            perm[i] = perm[belongTo] + nB[belongTo];
            nB[belongTo]++;
        }
    }
    t2 = sys_timer();           /* end angle method timer */
    *t_angle = t2 - t1;
    *pperm = perm;

    cleanCS(at);
    free(nB);
    free(jbuf);
    free(compress);

    return 0;
}