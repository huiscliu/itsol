
#include "solver-bicgstab.h"

#define  epsmac  1.0e-16

static double vec_norm_host(double *x, int n)
{
    int i;
    double t = 0.;
    
    assert(n >= 0);
    if (n > 0) assert(x != NULL);

    for (i = 0; i < n; i++)  t += x[i] * x[i];

    return sqrt(t);
}

static double vec_dot_host(double *x, double *y, int n)
{
    int i;
    double t = 0.;
    
    assert(n >= 0);
    if (n > 0) assert(x != NULL && y != NULL);

    for (i = 0; i < n; i++)  t += x[i] * y[i];

    return t;
}

int itsol_solver_bicgstab(ITS_SMat *Amat, ITS_PC *lu, double *bg, double *xg, double tol,
        int maxits, int *nits, double *res, FILE * fits)
{
    double *rg;
    double *rh;
    double *pg;
    double *ph;
    double *sg;
    double *sh;
    double *tg;
    double *vg;
    double *tp;
    double rho0 = 0, rho1 = 0;
    double alpha = 0, beta = 0, omega = 0;
    double residual;
    double err_rel = 0;
    int itr_out;
    int i;
    int n, retval = 0;

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

    residual = err_rel = vec_norm_host(rg, n);

    tol = residual * fabs(tol);

    for (itr_out = 0; itr_out < maxits; itr_out++) {
        rho1 = vec_dot_host(rg, rh, n);

        if (rho1 == 0) {
            fprintf(fits, "solver bicgstab failed.\n");
            break;
        }

        if (itr_out == 0) {

            for (i = 0; i < n; i++)
                pg[i] = rg[i];
        }
        else {
            beta = (rho1 * alpha) / (rho0 * omega);
            for (i = 0; i < n; i++) {
                pg[i] = rg[i] + beta * (pg[i] - omega * vg[i]);
            }
        }

        rho0 = rho1;

        /*  pc */
        if (lu == NULL) {
            memcpy(ph, pg, n * sizeof(double));
        }
        else {
            lu->precon(pg, ph, lu);
        }

        Amat->matvec(Amat, ph, vg);

        alpha = rho1 / vec_dot_host(rh, vg, n);
        for (i = 0; i < n; i++) {
            sg[i] = rg[i] - alpha * vg[i];
        }

        if (vec_norm_host(sg, n) <= 1e-60) {
            for (i = 0; i < n; i++) {
                xg[i] = xg[i] + alpha * ph[i];
            }

            Amat->matvec(Amat, xg, tp);
            for (i = 0; i < n; i++) rg[i] = bg[i] - tp[i];
            residual = vec_norm_host(rg, n);

            break;
        }

        if (lu == NULL) {
            memcpy(sh, sg, n * sizeof(double));
        }
        else {
            lu->precon(sg, sh, lu);
        }

        Amat->matvec(Amat, sh, tg);

        omega = vec_dot_host(tg, sg, n) / vec_dot_host(tg, tg, n);
        for (i = 0; i < n; i++) {
            xg[i] = xg[i] + alpha * ph[i] + omega * sh[i];
            rg[i] = sg[i] - omega * tg[i];
        }

        residual = vec_norm_host(rg, n);

        fprintf(fits, "it: %5d, abs res: %.6e, rel res: %.6e\n", itr_out, residual,
             (err_rel == 0 ? 0 : residual / err_rel));

        if (residual <= tol) break;
    }

    if (itr_out < maxits) itr_out += 1;

    free(rg);
    free(rh);
    free(pg);
    free(ph);
    free(sg);
    free(sh);
    free(tg);
    free(tp);
    free(vg);

    if (itr_out >= maxits) retval = 1;
    if (nits != NULL) *nits = itr_out;
    if (res != NULL) *res = residual;

    return retval;
}
