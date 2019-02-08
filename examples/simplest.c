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

    double *x = NULL, *rhs = NULL;
    int n, i;
    ITS_CooMat A;
    ITS_SOLVER s;

    /* case: COO formats */
    A = itsol_read_coo("pores3.coo");
    n = A.n;

    /* solution vectors */
    x = (double *)itsol_malloc(n * sizeof(double), "main");
    rhs = (double *)itsol_malloc(n * sizeof(double), "main");

    for( i = 0; i < n; i++ ) {
        x[i] = 0.0;
        rhs[i] = i;
    }

    /* create */
    itsol_solver_initialize(&s, ITS_PC_ARMS, &A);

    /* assemble */
    itsol_solver_assemble(&s);

    /* call itsol_solver_fgmres */
    itsol_solver_solve(&s, x, rhs);

    /* get results */
    printf("solver converged in %d steps...\n\n", s.nits);

    itsol_solver_finalize(&s);
    itsol_cleanCOO(&A);

    free(x);
    free(rhs);

    return 0;
}
