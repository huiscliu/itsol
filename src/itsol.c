
#include "itsol.h"

void itsol_solver_initialize(ITS_SOLVER *s, ITS_PC_TYPE pctype, ITS_CooMat *A)
{
    assert(s != NULL);
    assert(A != NULL);

    /* init */
    bzero(s, sizeof(*s));

    s->A = A;
    s->log = stdout;
    
    /* pc */
    s->pc_type = pctype;
    s->pc.log = s->log;
    itsol_pc_initialize(&s->pc, pctype);

    /* init parameters */
    itsol_solver_init_pars(&s->pars);
}

void itsol_solver_finalize(ITS_SOLVER *s)
{
    if (s == NULL) return;

    /* cleanup */
    if (s->csmat != NULL) itsol_cleanCS(s->csmat);
    s->csmat = NULL;

    itsol_pc_finalize(&s->pc);

    bzero(s, sizeof(*s));
}

int itsol_solver_assemble(ITS_SOLVER *s)
{
    ITS_PC_TYPE pctype;
    ITS_CooMat A;
    int ierr;
    FILE *log;

    assert(s != NULL);

    if (s->assembled) return 0;

    /* log */
    if (s->log == NULL) {
        log = stdout;
    }
    else {
        log = s->log;
    }

    /* assemble */
    pctype = s->pc_type;

    s->csmat = (ITS_SparMat *) itsol_malloc(sizeof(ITS_SparMat), "solver assemble");
    A = *s->A;

    if (pctype == ITS_PC_ILUC) {
        if ((ierr = itsol_COOcs(A.n, A.nnz, A.ma, A.ia, A.ja, s->csmat)) != 0) {
            fprintf(log, "solver assemble, COOcs error\n");
            return ierr;
        }

        /* smat */
        s->smat.n = A.n;
        s->smat.CS = s->csmat;               /* in column format */
        s->smat.matvec = itsol_matvecCSC;    /* column matvec */
    }
    else if(pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT || pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT
            || pctype == ITS_PC_ARMS) {
        if ((ierr = itsol_COOcs(A.n, A.nnz, A.ma, A.ja, A.ia, s->csmat)) != 0) {
            fprintf(log, "mainARMS: COOcs error\n");
            return ierr;
        }

        /* smat */
        s->smat.n = A.n;
        s->smat.CS = s->csmat;               /* in row format */
        s->smat.matvec = itsol_matvecCSR;    /* row matvec */
    }
    else {
        fprintf(log, "solver assemble, wrong preconditioner type\n");
        exit(-1);
    }

    itsol_pc_assemble(&s->pc);

    s->assembled = 1;
    return 0;
}

int itsol_solver_solve(ITS_SOLVER *s, double *x, double *rhs)
{
    FILE *log;
    ITS_PARS io;

    assert(s != NULL);
    assert(x != NULL);
    assert(rhs != NULL);

    io = s->pars;

    if (s->log == NULL) {
        log = stdout;
    }
    else {
        log = s->log;
    }

    return itsol_solver_fgmres(&s->smat, &s->pc, rhs, x, io.tol, io.restart, io.maxits, &s->nits, log);
}

void itsol_pc_initialize(ITS_PC *pc, ITS_PC_TYPE pctype)
{
    assert(pc != NULL);

    if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT) {
        pc->ILU = (ITS_ILUSpar *) itsol_malloc(sizeof(ITS_ILUSpar), "pc init");
    }
    else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
        pc->VBILU = (ITS_VBILUSpar *) itsol_malloc(sizeof(ITS_VBILUSpar), "pc init");
    }
    else if (pctype == ITS_PC_ARMS) {
        pc->ARMS = (ITS_ARMSpar *) itsol_malloc(sizeof(ITS_ARMSpar), "pc init");
    }
    else {
        fprintf(pc->log, "wrong preconditioner type\n");
        exit(-1);
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
