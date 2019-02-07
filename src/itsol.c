
#include "itsol.h"

void itsol_solver_initialize(ITS_SOLVER *s, ITS_PC_TYPE pctype, ITS_SparMat *A);
void itsol_solver_finalize(ITS_SOLVER *s);

void itsol_solver_assemble(ITS_SOLVER *s);

void itsol_pc_initialize(ITS_PC *pc, ITS_PC_TYPE pctype);
void itsol_pc_finalize(ITS_PC *pc);
void itsol_pc_assemble(ITS_PC *pc);
