/*-------------------------------------------------------------------*
 * main test driver for ILUK                                         *
 *-------------------------------------------------------------------*
 * Na Li, Oct 31, 2001                                               *
 *                                                                   *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu   *
 *-------------------------------------------------------------------*/
#include "itsol.h"

int main(void)
{
    int ierr = 0;
    /* --- plotting is for writing stats into ... in raw form */
    int plotting = 0, skip_its = 0;
    char pltfile[256];
    FILE *fits = NULL;
    /*-------------------- main structs and wraper structs.   */
    ITS_CsPtr csmat = NULL;         /* matrix in csr formt             */
    SMatptr MAT;                /* Matrix structure for matvecs    */
    SPreptr PRE;                /* general precond structure       */
    iluptr lu = NULL;           /* ilu preconditioner structure    */
    double *sol = NULL, *x = NULL, *rhs = NULL;
    /*-------------------- temp Harwell Boeing arrays */
    double *AA;
    int *JA, *IA;
    int n, nnz, lfil, rsa;
    /*-------------------- IO */
    FILE *flog = stdout, *fmat = NULL;
    io_t io;
    double tm1, tm2;
    int mat, numat, iparam, i;
    double terr;
    char line[ITS_MAX_LINE];

    MAT = (SMatptr) itsol_malloc(sizeof(SMat), "main:MAT");
    PRE = (SPreptr) itsol_malloc(sizeof(SPre), "main:PRE");

    /*------------------ read and set parameters and other inputs  */
    memset(&io, 0, sizeof(io));

    if (itsol_read_inputs("inputs", &io) != 0) {
        fprintf(flog, "Invalid inputs file...\n");
        exit(1);
    }

    /*------------------ file "matfile" contains paths to matrices */
    if (NULL == (fmat = fopen("matfile", "r"))) {
        fprintf(flog, "Can't open matfile...\n");
        exit(2);
    }

    memset(line, 0, ITS_MAX_LINE);
    fgets(line, ITS_MAX_LINE, fmat);

    if ((numat = atoi(line)) <= 0) {
        fprintf(flog, "Invalid count of matrices...\n");
        exit(3);
    }

    /*-------------------- open file ILUK.out for all performance
      results of this run (all matrices and params) 
      also set io->PrecMeth */
    strcpy(io.outfile, "ILUK.out");
    strcpy(io.PrecMeth, "ILUK");

    if (NULL == (io.fout = fopen(io.outfile, "w"))) {
        fprintf(flog, "Can't open output file %s...\n", io.outfile);
        exit(4);
    }

    /* ------------------- LOOP THROUGH MATRICES ------------------- */
    for (mat = 1; mat <= numat; mat++) {
        if (itsol_get_matrix_info(fmat, &io) != 0) {
            fprintf(flog, "Invalid format in matfile_hb...\n");
            exit(5);
        }

        fprintf(flog, "MATRIX: %s...\n", io.MatNam);

        /*------------------------- Read matrix */
        csmat = (ITS_CsPtr) itsol_malloc(sizeof(SparMat), "main");

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

        /*-------------------- free COO/HB arrays  */
        free(AA);
        AA = NULL;
        free(JA);
        JA = NULL;
        free(IA);
        IA = NULL;

        /*---------------------------------------------------------*/
        x = (double *)itsol_malloc(io.ndim * sizeof(double), "main");
        itsol_output_header(&io);

        /*-------------------- set initial lfil and tol */
        lfil = io.fill_lev;
        io.tol0 = 0.0;

        /*--------------------------LOOP THROUGH PARAMETERS */
        for (iparam = 1; iparam <= io.nparam; iparam++) {
            fprintf(flog, "iparam = %d\n", iparam);
            lu = (iluptr) itsol_malloc(sizeof(ILUSpar), "main");
            fprintf(flog, "begin iluk(%d)\n", lfil);
            tm1 = itsol_get_time();

            /*-------------------- call ILUK preconditioner set-up  */
            ierr = itsol_pc_ilukC(lfil, csmat, lu, flog);

            /*----------------------------------------------------- */
            tm2 = itsol_get_time();

            if (ierr == -2) {
                fprintf(io.fout, "zero diagonal element found...\n");
                itsol_cleanILU(lu);
                goto NEXT_MAT;
            }
            else if (ierr != 0) {
                fprintf(flog, "*** iluk error, ierr != 0 ***\n");
                exit(-1);
            }

            io.tm_p = tm2 - tm1;
            io.fillfact = itsol_nnz_ilu(lu) / (double)(io.nnz + 1);
            fprintf(flog, "iluk ends, fill factor (mem used) = %f\n", io.fillfact);

            if (skip_its) {
                io.its = -1;
                io.tm_i = -1;
                io.enorm = -1;
                io.rnorm = -1;
                goto NEXT_PARA;
            }

            /*------------------------ condition number estimation */
            if (itsol_condestLU(lu, flog) != 0) {
                fprintf(flog, "Not attempting iterative solution.\n");
                fprintf(io.fout, "Not attempting iterative solution.\n");
                io.its = -1;
                io.tm_i = -1;
                io.enorm = -1;
                io.rnorm = -1;
                goto NEXT_PARA;
            }

            /* ----------------------- initial guess */
            /*for( i = 0; i < io.ndim; i++ ) x[i] = 0.0; */
            itsol_randvec(x, n);

            /*-------------------- create a file for printing
              'its -- time -- res' info from fgmres */
            if (plotting) {
                sprintf(pltfile, "%s_ILUK_F%05d", io.MatNam, lfil);
                if (NULL == (fits = fopen(pltfile, "w"))) {
                    fprintf(flog, "Can't open output file %s...\n", pltfile);
                    exit(8);
                }
            }
            else
                fits = NULL;

            /*-------------------- set up the structs before calling itsol_solver_fgmres */
            MAT->n = n;
            MAT->CS = csmat;
            MAT->matvec = itsol_matvecCSR;
            PRE->ILU = lu;
            PRE->precon = itsol_preconILU;

            /*-------------------- call itsol_solver_fgmres */
            io.its = io.maxits;
            tm1 = itsol_get_time();
            itsol_solver_fgmres(MAT, PRE, rhs, x, io.tol, io.im, &io.its, fits);
            tm2 = itsol_get_time();
            io.tm_i = tm2 - tm1;

            if (io.its < io.maxits)
                fprintf(flog, "param %03d OK: converged in %d steps...\n\n", iparam, io.its);
            else
                fprintf(flog, "not converged in %d steps...\n\n", io.maxits);

            if (fits) fclose(fits);

            /*-------------------- calculate error norm */
            terr = 0.0;
            for (i = 0; i < io.ndim; i++) terr += (sol[i] - x[i]) * (sol[i] - x[i]);
            io.enorm = sqrt(terr);

            /*-------------------- calculate residual norm */
            itsol_matvec(csmat, x, sol);
            terr = 0.0;

            for (i = 0; i < io.ndim; i++) terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);

            io.rnorm = sqrt(terr);

            /*-------------------- Test with next params   */
NEXT_PARA:
            itsol_output_result(lfil, &io, iparam);
            lfil += io.fill_lev_inc;
            itsol_cleanILU(lu);
        }

        /*-------------------- Test with next matrix   */
NEXT_MAT:
        itsol_cleanCS(csmat);
        free(sol);
        free(x);
        free(rhs);
    }

    fclose(io.fout);
    if (flog != stdout) fclose(flog);
    fclose(fmat);
    free(MAT);
    free(PRE);
    return 0;
}
