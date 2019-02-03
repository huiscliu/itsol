/*-----------------------------------------------------------------*
 * main test driver for VBILUT                                     *
 *-----------------------------------------------------------------*
 * Na Li, Aug 26, 2001 -- YS 2005                                  *
 *                                                                 *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu *
 *-----------------------------------------------------------------*/
#include "itsol.h"

int main(void)
{
    int ierr = 0;
    /*-------------------------------------------------------------------
     * options
     *-----------------------------------------------------------------*/
    int plotting = 0, output_mat = 0, skip_its = 0, diagscal = 0;
    char pltfile[256];
    FILE *fits = NULL;
    /*-------------------- main structs and wraper structs.   */
    csptr csmat = NULL;         /* matrix in csr formt             */
    vbsptr vbmat = NULL;
    vbiluptr lu = NULL;         /* vbilu preconditioner structure  */
    SMatptr MAT;                /* Matrix structure for matvecs    */
    SPreptr PRE;                /* general precond structure       */
    double *sol = NULL, *x = NULL, *prhs = NULL, *rhs = NULL;
    /*-------------------- temp COO/HB arrays - fortran style */
    int *IA, *JA;
    double *AA;
    /*---------------------------------------------------------*/
    int n, nnz, rsa;
    BData *w = NULL;
    int lfil, max_blk_sz = MAX_BLOCK_SIZE * MAX_BLOCK_SIZE * sizeof(double);
    int nBlock, *nB = NULL, *perm = NULL;
    double tol;
    /*-------------------- IO */
    FILE *flog = stdout, *fmat = NULL;
    io_t io;
    /*---------------------------------------------------------*/
    double tm1, tm2;
    int mat, numat, iparam, i;
    double terr;
    char line[MAX_LINE];

    MAT = (SMatptr) itsol_malloc(sizeof(SMat), "main:MAT");
    PRE = (SPreptr) itsol_malloc(sizeof(SPre), "main:PRE");

    /*------------------ read and set parameters and other inputs  */
    memset(&io, 0, sizeof(io));

    if (itsol_read_inputs("inputs", &io) != 0) {
        fprintf(flog, "Invalid inputs file...\n");
        exit(1);
    }

    /*------------------ set any parameters manually */
    /* io.eps  is the angle tolerance for grouping two columns in same
       supernode. This is a cosine and should be <= 1.  */
    io.eps = 0.8;

    /*------------------ file "matfile" contains paths to matrices */
    if (NULL == (fmat = fopen("matfile", "r"))) {
        fprintf(flog, "Can't open matfile...\n");
        exit(2);
    }

    memset(line, 0, MAX_LINE);
    fgets(line, MAX_LINE, fmat);

    if ((numat = atoi(line)) <= 0) {
        fprintf(flog, "Invalid count of matrices...\n");
        exit(3);
    }

    /*-------------------- open file ILUT.out for all performance
      results of this run (all matrices and params) 
      also set io->PrecMeth */
    strcpy(io.outfile, "VBILUT.out");
    strcpy(io.PrecMeth, "Variable Block ILUT (VBILUT)");

    if (NULL == (io.fout = fopen(io.outfile, "w"))) {
        fprintf(flog, "Can't open output file %s...\n", io.outfile);
        exit(4);
    }

    /*-------------------- LOOP THROUGH MATRICES --------------------*/
    for (mat = 1; mat <= numat; mat++) {
        if (itsol_get_matrix_info(fmat, &io) != 0) {
            fprintf(flog, "Invalid format in matfile_hb...\n");
            exit(5);
        }

        fprintf(flog, "MATRIX: %s...\n", io.MatNam);

        /* ------------------- Read in matrix and allocate memory-------- */
        csmat = (csptr) itsol_malloc(sizeof(SparMat), "main");

        /*-------------------- case: COO formats */
        if (io.Fmt > HB) {
            ierr = itsol_read_coo(&AA, &JA, &IA, &io, &rhs, &sol, 0);
            if (ierr == 0)
                fprintf(flog, "matrix read successfully\n");
            else {
                fprintf(flog, "read_coo error = %d\n", ierr);
                exit(6);
            }

            n = io.ndim;
            nnz = io.nnz;

            /*-------------------- conversion from COO to CSR format */
            if ((ierr = itsol_COOcs(n, nnz, AA, JA, IA, csmat)) != 0) {
                fprintf(stderr, "mainARMS: COOcs error\n");
                return ierr;
            }
        }
        else if (io.Fmt == HB) {
            /*-------------------- NOTE: (AA,JA,IA) is in CSR format */
            ierr = itsol_readhb_c(&n, &AA, &JA, &IA, &io, &rhs, &sol, &rsa);

            if (ierr != 0) {
                fprintf(flog, "readhb_c error = %d\n", ierr);
                exit(7);
            }

            nnz = io.nnz;

            if ((ierr = itsol_CSRcs(n, AA, JA, IA, csmat, rsa)) != 0) {
                fprintf(flog, "readhb_c: CSRcs error\n");
                return ierr;
            }
        }

        /*----------------------- Free Memory */
        free(AA);
        AA = NULL;
        free(JA);
        JA = NULL;
        free(IA);
        IA = NULL;

        /*----------------------- Daigonal Scaling */
        if (diagscal == 1) {
            int nrm = 1;
            double *diag;

            diag = (double *)itsol_malloc(sizeof(double) * n, "mainILUC:diag");
            ierr = itsol_roscalC(csmat, diag, nrm);

            if (ierr != 0) {
                fprintf(stderr, "main-vbiluk: roscal: a zero row...\n");
                return ierr;
            }

            ierr = itsol_coscalC(csmat, diag, nrm);

            if (ierr != 0) {
                fprintf(stderr, "main-vbiluk: roscal: a zero col...\n");
                return ierr;
            }

            free(diag);
        }

        /*------------------------------------------------------------*/
        x = (double *)itsol_malloc(io.ndim * sizeof(double), "main");
        ierr = itsol_init_blocks(csmat, &nBlock, &nB, &perm, io.eps, &io.tm_h, &io.tm_a);
        io.tm_b = io.tm_h + io.tm_a;

        if (ierr != 0) {
            fprintf(flog, "*** in init_blocks ierr != 0 ***\n");
            exit(8);
        }

        /*------------- permutes the rows and columns of the matrix */
        if (itsol_dpermC(csmat, perm) != 0) {
            fprintf(flog, "*** dpermC error ***\n");
            exit(9);
        }

        /*------------- permutes right hand side */
        prhs = (double *)itsol_malloc(n * sizeof(double), "main");

        for (i = 0; i < n; i++) prhs[perm[i]] = rhs[i];

        /*-------------------- convert to block matrix. */
        vbmat = (vbsptr) itsol_malloc(sizeof(VBSparMat), "main");
        ierr = itsol_csrvbsrC(1, nBlock, nB, csmat, vbmat);

        if (ierr != 0) {
            fprintf(flog, "*** in csrvbsr ierr != 0 ***\n");
            exit(10);
        }

        /*------------------------- OUTPUT MATRIX */
        if (output_mat) {
            char matdata[MAX_LINE];
            FILE *fmatlab;
            int ii, jj;

            sprintf(matdata, "%s.dat", io.MatNam);
            if (NULL != (fmatlab = fopen(matdata, "w"))) {
                fprintf(fmatlab, "%d %d 0\n", csmat->n, csmat->n);
                for (ii = 0; ii < csmat->n; ii++)
                    for (jj = 0; jj < csmat->nzcount[ii]; jj++)
                        fprintf(fmatlab, "%d %d 1\n", ii + 1, csmat->ja[ii][jj] + 1);
                fclose(fmatlab);
            }
        }

        /*------------- info of Block matrix */
        io.rt_v = (double)csmat->n / (double)vbmat->n;
        io.rt_e = itsol_nnz_cs(csmat) / 1. / itsol_nnzVBMat(vbmat);
        io.ceff = itsol_nnz_cs(csmat) / 1. / itsol_memVBMat(vbmat) * 100;

        /*---------------------------*/
        itsol_output_header_vb(&io);
        lfil = io.lfil0;
        tol = io.tol0;
        w = (BData *) itsol_malloc(vbmat->n * sizeof(BData), "main");

        for (i = 0; i < vbmat->n; i++)
            w[i] = (double *)itsol_malloc(max_blk_sz, "main");

        /* ----------------------- LOOP THROUGH PARAMETERS ------------- */
        for (iparam = 1; iparam <= io.nparam; iparam++) {
            fprintf(flog, "iparam = %d\n", iparam);
            lu = (vbiluptr) itsol_malloc(sizeof(VBILUSpar), "main");
            fprintf(flog, "begin vbilut\n");

            tm1 = itsol_get_time();

            /*-------------------- call VBILUT preconditioner set-up  */
            ierr = itsol_pc_vbilutC(vbmat, lu, lfil, tol, w, flog);

            /*----------------------------------------------------- */
            tm2 = itsol_get_time();

            if (ierr == -2) {
                fprintf(io.fout, "Singular diagonal block...\n");
                itsol_cleanVBILU(lu);
                goto NEXT_MAT;
            }
            else if (ierr != 0) {
                fprintf(flog, "*** vbilu error, ierr != 0 ***\n");
                exit(11);
            }

            io.tm_p = tm2 - tm1;
            io.fillfact = itsol_nnz_vbilu(lu) / (double)(io.nnz + 1);
            fprintf(flog, "vbilut ends, fill factor (mem used) = %f\n", io.fillfact);

            /*----------------------------------------------------- */
            if (skip_its) {
                io.its = -1;
                io.tm_i = -1;
                io.enorm = -1;
                io.rnorm = -1;
                goto NEXT_PARA;
            }

            /*------------- get rough idea of cond number - exit if too big */
            if (itsol_VBcondestC(lu, flog) != 0) {
                fprintf(flog, "Not attempting iterative solution.\n");
                fprintf(io.fout, "Not attempting iterative solution.\n");
                io.its = -1;
                io.tm_i = -1;
                io.enorm = -1;
                io.rnorm = -1;
                goto NEXT_PARA;
            }

            /*-------------------- initial guess */
            /* for( i = 0; i < io.ndim; i++ ) x[i] = 0.0; */
            itsol_randvec(x, n);

            /*-------------------- create a file for printing
              'its -- time -- res' info from fgmres */
            if (plotting) {
                sprintf(pltfile, "%s_VBILUT_F%05d_T%08.6f", io.MatNam, lfil, tol);
                if (NULL == (fits = fopen(pltfile, "w"))) {
                    fprintf(flog, "Can't open output file %s...\n", pltfile);
                    exit(12);
                }
            }
            else
                fits = NULL;

            /*-------------------- set up the structs before calling itsol_solver_fgmres */
            MAT->n = n;
            MAT->CS = csmat;
            MAT->matvec = itsol_matvecCSR;
            PRE->VBILU = lu;
            PRE->precon = itsol_preconVBR;

            /*-------------------- call itsol_solver_fgmres */
            io.its = io.maxits;
            tm1 = itsol_get_time();
            itsol_solver_fgmres(MAT, PRE, prhs, x, io.tol, io.im, &io.its, fits);
            tm2 = itsol_get_time();
            io.tm_i = tm2 - tm1;

            if (io.its < io.maxits)
                fprintf(flog, "param %03d OK: converged in %d steps...\n\n", iparam, io.its);
            else
                fprintf(flog, "not converged in %d steps...\n\n", io.maxits);

            if (fits) fclose(fits);

            /*---------------------- calculate error norm */
            /*---------------------- P*sol ?= x */
            terr = 0.0;
            for (i = 0; i < io.ndim; i++)
                terr += (sol[i] - x[perm[i]]) * (sol[i] - x[perm[i]]);

            io.enorm = sqrt(terr);

            /*---------------------- calculate residual norm */
            itsol_vbmatvec(vbmat, x, sol);

            terr = 0.0;
            for (i = 0; i < io.ndim; i++)
                terr += (prhs[i] - sol[i]) * (prhs[i] - sol[i]);

            io.rnorm = sqrt(terr);

            /*-------------------- next params */
NEXT_PARA:
            itsol_output_result(lfil, &io, iparam);
            lfil += io.lfilInc;
            tol *= io.tolMul;
            itsol_cleanVBILU(lu);
        }

NEXT_MAT:
        /*  output_blocks( nBlock, nB, io.fout ); */
        for (i = 0; i < vbmat->n; i++) free(w[i]);

        free(w);
        itsol_cleanCS(csmat);
        itsol_cleanVBMat(vbmat);
        free(nB);
        free(perm);
        free(sol);
        free(x);
        free(prhs);
        free(rhs);
    }

    fclose(io.fout);
    if (flog != stdout) fclose(flog);
    fclose(fmat);

    free(MAT);
    free(PRE);

    return 0;
}
