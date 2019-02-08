
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

    /* pc assemble */
    itsol_pc_assemble(s);

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

    pc->pc_type = pctype;

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
    ITS_PC_TYPE pctype;

    if (pc == NULL) return;

    pctype = pc->pc_type;
    if (pctype == ITS_PC_ILUC || pctype == ITS_PC_ILUK || pctype == ITS_PC_ILUT) {
        itsol_cleanILU(pc->ILU);
        pc->ILU = NULL;
    }
    else if (pctype == ITS_PC_VBILUK || pctype == ITS_PC_VBILUT) {
        itsol_cleanVBILU(pc->VBILU);
        pc->VBILU = NULL;
    }
    else if (pctype == ITS_PC_ARMS) {
        itsol_cleanARMS(pc->ARMS);
        pc->ARMS = NULL;
    }
    else {
        fprintf(pc->log, "wrong preconditioner type\n");
        exit(-1);
    }
}

int itsol_pc_assemble(ITS_SOLVER *s)
{
    ITS_PC_TYPE pctype;
    int ierr;
    ITS_PARS p;
    ITS_PC *pc;

    assert(s != NULL);
    pc = &s->pc;

    /* type */
    pctype = pc->pc_type;
    p = s->pars;

    if (pctype == ITS_PC_ILUC) {
        pc->precon = itsol_preconLDU;
    }
    else if (pctype == ITS_PC_ILUK) {
        ierr = itsol_pc_ilukC(p.fill_lev, s->csmat, pc->ILU, pc->log);

        if (ierr != 0) {
            fprintf(pc->log, "pc assemble, ILUK error\n");
            return ierr;
        }

        pc->precon = itsol_preconILU;
    }
    else if (pctype == ITS_PC_ILUT) {
        ierr = itsol_pc_ilut(s->csmat, pc->ILU, p.lfil0, p.tol0, pc->log);

        if (ierr != 0) {
            fprintf(pc->log, "pc assemble, ILUK error\n");
            return ierr;
        }

        pc->precon = itsol_preconILU;
    }
    else if (pctype == ITS_PC_VBILUK) {
        pc->precon = itsol_preconVBR;
    }
    else if (pctype == ITS_PC_VBILUT) {
        pc->precon = itsol_preconVBR;
    }
    else if (pctype == ITS_PC_ARMS) {
        pc->precon = itsol_preconARMS;
    }
    else {
        fprintf(pc->log, "wrong preconditioner type\n");
        exit(-1);
    }

    return 0;
}


void itsol_solver_set_pars(ITS_SOLVER *s, ITS_PARS par)
{
    ITS_PARS *p;

    assert(s != NULL);

    memcpy(&s->pars, &par, sizeof(par));

    /* update arms pars */
    p = &s->pars;
    itsol_set_arms_pars(p, p->diagscal, p->ipar, p->dropcoef, p->lfil_arr);
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

    /* arms */
    p->diagscal = 1;

    /* init arms pars */
    itsol_set_arms_pars(p, p->diagscal, p->ipar, p->dropcoef, p->lfil_arr);
}
