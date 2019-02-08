
#include "itsol.h"

void itsol_solver_initialize(ITS_SOLVER *s, ITS_PC_TYPE pctype, ITS_CooMat *A)
{
    assert(s != NULL);
    assert(A != NULL);

    /* init */
    bzero(s, sizeof(*s));

    s->A = A;
    
    /* pc */
    s->pc_type = pctype;
    itsol_pc_initialize(&s->pc, pctype);

    /* init parameters */
    itsol_solver_init_pars(&s->pars);
}

void itsol_solver_finalize(ITS_SOLVER *s)
{
    if (s == NULL) return;

    itsol_pc_finalize(&s->pc);
}

void itsol_solver_assemble(ITS_SOLVER *s)
{
    assert(s != NULL);

    itsol_pc_assemble(&s->pc);
}

void itsol_solver_solve(ITS_SOLVER *s, double *x, double *rhs);

void itsol_pc_initialize(ITS_PC *pc, ITS_PC_TYPE pctype)
{
    assert(pc != NULL);

    if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT) {
    }
    else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
    }
    else if (pctype == ITS_PC_ARMS) {
    }
    else {
    }
}

void itsol_pc_finalize(ITS_PC *pc)
{
    if (pc == NULL) return;
}

void itsol_pc_assemble(ITS_PC *pc)
{
    if (pc == NULL) return;
}


void itsol_solver_set_pars(ITS_SOLVER *s, ITS_PARS par)
{
    assert(s != NULL);

    memcpy(&s->pars, &par, sizeof(par));
}

void itsol_solver_init_pars(ITS_PARS *p)
{
    assert(p != NULL);

    /* parameters from inputs -----------------------------------------*/
    p->restart = 30;               /* Dim of Krylov subspace [fgmr]   */
    p->maxits = 1000;              /* maximum number of fgmres iters  */
    p->tol = 1e-6;                 /* tolerance for stopping fgmres   */

    p->eps = 0.8;
    p->lfil0 = 50;                 /* initial lfil                    */
    p->tol0 = 1e-3;                /* initial drop tolerance          */
    p->fill_lev = 1;               /* initial level of fill for ILUK  */

    /* value always set to 1           */
    p->perm_type = 0;              /* indset perms (0) or PQ perms (1)*/
    p->Bsize = 30;                 /* block size - dual role. see input file */
}
