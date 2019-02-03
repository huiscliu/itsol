
#ifndef __ITSOL_ILUK_H__
#define __ITSOL_ILUK_H__

#include "type-defs.h"
#include "sets.h"

#ifdef __cplusplus
extern "C" {
#endif

int lofC(int lofM, csptr csmat, iluptr lu, FILE *fp); 
int ilukC(int lofM, csptr csmat, iluptr lu, FILE *fp);

#ifdef __cplusplus
}
#endif
#endif
