
#include "solver-bicgstab.h"

#define  epsmac  1.0e-16

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

    mv_amxpbyz_host(A, xg, -1, bg, 1, rg);
    for (i = 0; i < n; i++) {
        rh[i] = rg[i];
        sh[i] = ph[i] = 0.;
    }

    residual = err_rel = vec_norm_host(rg, n);

    tol = residual * fabs(tol);

    for (itr_out = 0; itr_out < maxits; itr_out++) {
        rho1 = vec_dot_host(rg, rh, n);

        if (rho1 == 0) {
            solver_lprintf("BICTSTAB, CPU: method failed.!\n");
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

        if (pc.type != PC_NON) {
            vec_set_value_host(ph, 0., n);
            pc.solve(pc, ph, pg);
        }
        else {
            for (i = 0; i < n; i++) {
                ph[i] = pg[i];
            }
        }

        mv_amxpbyz_host(A, ph, 1, pg, 0, vg);

        alpha = rho1 / vec_dot_host(rh, vg, n);
        for (i = 0; i < n; i++) {
            sg[i] = rg[i] - alpha * vg[i];
        }

        if (vec_norm_host(sg, n) <= PASS_SOLVER_BREAKDOWN) {
            solver_lprintf
                ("BICTSTAB, CPU: ||s|| is too small: %f, terminated.\n",
                 vec_norm_host(sg, n));

            for (i = 0; i < n; i++) {
                xg[i] = xg[i] + alpha * ph[i];
            }

            mv_amxpbyz_host(A, xg, -1, bg, 1, rg);
            residual = vec_norm_host(rg, n);

            break;
        }

        if (pc.type != PC_NON) {
            vec_set_value_host(sh, 0., n);
            pc.solve(pc, sh, sg);
        }

        mv_amxpbyz_host(A, sh, 1, pg, 0, tg);

        omega = vec_dot_host(tg, sg, n) / vec_dot_host(tg, tg, n);
        for (i = 0; i < n; i++) {
            xg[i] = xg[i] + alpha * ph[i] + omega * sh[i];
            rg[i] = sg[i] - omega * tg[i];
        }

        residual = vec_norm_host(rg, n);

        frintf(fits, "it: %5d, abs res: %.6e, rel res: %.6e\n", itr_out, residual,
             (err_rel == 0 ? 0 : residual / err_rel));

        if (residual <= tol) break;
    }

    if (itr_out < maxits) itr_out += 1;

end:

    solver.residual = residual;
    solver.nits = itr_out;

    array_free_host < double > (rg);
    array_free_host < double > (rh);
    array_free_host < double > (pg);
    array_free_host < double > (ph);
    array_free_host < double > (sg);
    array_free_host < double > (sh);
    array_free_host < double > (tg);
    array_free_host < double > (vg);
    array_free_host < double > (rscl);

    if (iter >= maxits) retval = 1;
    if (nits != NULL) *nits = itr_out;

    return retval;
}
