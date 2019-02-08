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

    double *sol = NULL, *x = NULL, *rhs = NULL;
    int n;
    int i;
    double terr, norm;
    ITS_CooMat A;
    int its;
    ITS_SparMat *csmat = NULL;
    ITS_SOLVER s;

    /*-------------------- case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;

    /*---------------------------------------------------------*/
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");
    sol = (double *)itsol_malloc(n * sizeof(double), "main");

    for( i = 0; i < n; i++ ) x[i] = 0.0;
    for (i = 0; i < n; i++) rhs[i] = i;

    /*-------------------- call itsol_solver_fgmres */
    itsol_solver_initialize(&s, ITS_PC_ILUK, &A);

    itsol_solver_assemble(&s);

    itsol_solver_solve(&s, x, rhs);

    its = s.nits;
    printf("solver converged in %d steps...\n\n", its);

    /*-------------------- calculate residual norm */
    csmat = s.csmat;
    itsol_matvec(csmat, x, sol);

    /* error */
    terr = 0.0;
    norm = 0.;
    for (i = 0; i < A.n; i++) {
        terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);

        norm += rhs[i] * rhs[i];
    }

    printf("residual: %e, relative residual: %e\n\n", sqrt(terr), sqrt(terr / norm));

    itsol_solver_finalize(&s);
    itsol_cleanCOO(&A);

    free(sol);
    free(x);
    free(rhs);

    return 0;
}
