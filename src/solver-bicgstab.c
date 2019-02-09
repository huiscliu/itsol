
#include "solver-bicgstab.h"

#define  epsmac  1.0e-16

int itsol_solver_bicgstab(ITS_SMat *Amat, ITS_PC *lu, double *rhs, double *sol, double tol,
        int maxits, int *nits, double *res, FILE * fits)
{
    int n = Amat->n;
    int retval, i = 0;

    retval = 0;

    for (i = 0; i < maxits; i++) {
        //Amat->matvec(Amat, sol, vv);
    }

    if (i >= maxits) retval = 1;
    *nits = i;

    return retval;
}
