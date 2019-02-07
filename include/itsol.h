
#ifndef __ITSOL_INCLUDED_PROTOS_H__
#define __ITSOL_INCLUDED_PROTOS_H__

#include "auxill.h"
#include "solver-fgmres.h"

#include "pc-arms2.h"
#include "pc-iluk.h"
#include "pc-ilutc.h"
#include "pc-ilut.h"
#include "pc-ilutpc.h"
#include "pc-pilu.h"
#include "pc-vbiluk.h"
#include "pc-vbilut.h"

#ifdef __cplusplus
extern "C" {
#endif

void itsol_solver_initialize(ITS_SOLVER *s, ITS_PC_TYPE pctype, ITS_SparMat *A);
void itsol_solver_finalize(ITS_SOLVER *s);

void itsol_pc_initialize(ITS_PC *pc, ITS_PC_TYPE pctype);
void itsol_pc_finalize(ITS_PC *pc);

double * itsol_solver_get_solution(ITS_SOLVER *s);

#ifdef __cplusplus
}
#endif
#endif 
