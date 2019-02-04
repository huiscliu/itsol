
#ifndef __ITSOL_ILUK_H__
#define __ITSOL_ILUK_H__

#include "type-defs.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

int itsol_pc_lofC(int lofM, ITS_CsPtr csmat, ITS_IluPtr lu, FILE *fp); 
int itsol_pc_ilukC(int lofM, ITS_CsPtr csmat, ITS_IluPtr lu, FILE *fp);

#ifdef __cplusplus
}
#endif
#endif
