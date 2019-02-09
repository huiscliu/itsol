
#include "solver-bicgstab.h"

int itsol_solver_bicgstab(ITS_SMat *Amat, ITS_PC *lu, double *bg, double *xg, double tol,
        int maxits, int *nits, double *res, FILE * fits)
{
    double *rg, *rh, *pg, *ph, *sg, *sh, *tg, *vg, *tp;
    double r0 = 0, r1 = 0, pra = 0, prb = 0, prc = 0;
    double residual, err_rel = 0;
    int i, n, retval = 0;
    int itr = 0.;

    n = Amat->n;
    rg = itsol_malloc(n * sizeof(double), "bicgstab");
    rh = itsol_malloc(n * sizeof(double), "bicgstab");
    pg = itsol_malloc(n * sizeof(double), "bicgstab");
    ph = itsol_malloc(n * sizeof(double), "bicgstab");
    sg = itsol_malloc(n * sizeof(double), "bicgstab");
    sh = itsol_malloc(n * sizeof(double), "bicgstab");
    tg = itsol_malloc(n * sizeof(double), "bicgstab");
    vg = itsol_malloc(n * sizeof(double), "bicgstab");
    tp = itsol_malloc(n * sizeof(double), "bicgstab");

    Amat->matvec(Amat, xg, tp);
    for (i = 0; i < n; i++) rg[i] = bg[i] - tp[i];

    for (i = 0; i < n; i++) {
        rh[i] = rg[i];
        sh[i] = ph[i] = 0.;
    }

    residual = err_rel = itsol_norm(rg, n);
    tol = residual * fabs(tol);

    if (tol == 0.) goto skip;

    for (itr = 0; itr < maxits; itr++) {
        r1 = itsol_dot(rg, rh, n);

        if (r1 == 0) {
            fprintf(fits, "solver bicgstab failed.\n");
            break;
        }

        if (itr == 0) {

            for (i = 0; i < n; i++)
                pg[i] = rg[i];
        }
        else {
            prb = (r1 * pra) / (r0 * prc);
            for (i = 0; i < n; i++) {
                pg[i] = rg[i] + prb * (pg[i] - prc * vg[i]);
            }
        }

        r0 = r1;

        /*  pc */
        if (lu == NULL) {
            memcpy(ph, pg, n * sizeof(double));
        }
        else {
            lu->precon(pg, ph, lu);
        }

        Amat->matvec(Amat, ph, vg);

        pra = r1 / itsol_dot(rh, vg, n);
        for (i = 0; i < n; i++) {
            sg[i] = rg[i] - pra * vg[i];
        }

        if (itsol_norm(sg, n) <= 1e-60) {
            for (i = 0; i < n; i++) {
                xg[i] = xg[i] + pra * ph[i];
            }

            Amat->matvec(Amat, xg, tp);
            for (i = 0; i < n; i++) rg[i] = bg[i] - tp[i];
            residual = itsol_norm(rg, n);

            break;
        }

        if (lu == NULL) {
            memcpy(sh, sg, n * sizeof(double));
        }
        else {
            lu->precon(sg, sh, lu);
        }

        Amat->matvec(Amat, sh, tg);

        prc = itsol_dot(tg, sg, n) / itsol_dot(tg, tg, n);
        for (i = 0; i < n; i++) {
            xg[i] = xg[i] + pra * ph[i] + prc * sh[i];
            rg[i] = sg[i] - prc * tg[i];
        }

        residual = itsol_norm(rg, n);

        fprintf(fits, "%8d   %10.2e\n", itr, residual / err_rel);

        if (residual <= tol) break;
    }

    if (itr < maxits) itr += 1;

skip:
    free(rg);
    free(rh);
    free(pg);
    free(ph);
    free(sg);
    free(sh);
    free(tg);
    free(tp);
    free(vg);

    if (itr >= maxits) retval = 1;
    if (nits != NULL) *nits = itr;
    if (res != NULL) *res = residual;

    return retval;
}
